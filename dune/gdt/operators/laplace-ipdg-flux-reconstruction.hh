// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_OPERATORS_LAPLACE_IPDG_FLUX_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_LAPLACE_IPDG_FLUX_RECONSTRUCTION_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/orthonormal.hh>
#include <dune/gdt/local/integrands/ipdg.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/print.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>

namespace Dune {
namespace GDT {


/**
 * \attention Handles all boundary intersections as Dirichlet boundary intersections!
 *
 * \todo Add volume Dofs for higher order analoggously to the RT interpolation!
 *
 * \todo Directly make use of LocalIPDGIntegrands::InnerPenalty and LocalLaplaceIPDGIntegrands::InnerCoupling, (not
 *       clear how yet)!
 */
template <class M, class AssemblyGridView, class SGV = AssemblyGridView, class RGV = AssemblyGridView>
class LaplaceIpdgFluxReconstructionOperator : public OperatorInterface<M, SGV, 1, 1, RGV::dimension, 1, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  using BaseType = OperatorInterface<M, SGV, 1, 1, RGV::dimension, 1, RGV>;
  using ThisType = LaplaceIpdgFluxReconstructionOperator;

public:
  using typename BaseType::F;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AGV = AssemblyGridView;
  using AssemblyGridViewType = AGV;
  using E = XT::Grid::extract_entity_t<AssemblyGridViewType>;
  using I = XT::Grid::extract_intersection_t<AssemblyGridViewType>;
  static constexpr size_t d = RGV::dimension;

  LaplaceIpdgFluxReconstructionOperator(AssemblyGridViewType assembly_grid_view,
                                        const SourceSpaceType& src_spc,
                                        const RangeSpaceType& rng_spc,
                                        const double& symmetry_prefactor,
                                        const double& inner_penalty,
                                        const double& dirichlet_penalty,
                                        XT::Functions::GridFunction<E, d, d> diffusion,
                                        XT::Functions::GridFunction<E, d, d> weight_function = {1.},
                                        const std::function<double(const I&)>& intersection_diameter =
                                            LocalIPDGIntegrands::internal::default_intersection_diameter<I>(),
                                        const std::string& logging_prefix = "")
    : BaseType({},
               logging_prefix.empty() ? "LaplaceIpdgFluxReconstructionOperator" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , assembly_grid_view_(assembly_grid_view)
    , source_space_(src_spc)
    , range_space_(rng_spc)
    , symmetry_prefactor_(symmetry_prefactor)
    , inner_penalty_(inner_penalty)
    , dirichlet_penalty_(dirichlet_penalty)
    , diffusion_(diffusion)
    , weight_function_(weight_function)
    , intersection_diameter_(intersection_diameter)
    , element_mapper_(assembly_grid_view_)
  {
    DUNE_THROW_IF(range_space_.type() != SpaceType::raviart_thomas, Exceptions::operator_error, "");
    DUNE_THROW_IF(range_space_.max_polorder() != 0, Exceptions::operator_error, "Not implemented yet!");
  }

  // manual copy ctor due to element_mapper_
  LaplaceIpdgFluxReconstructionOperator(const ThisType& other)
    : BaseType(other)
    , assembly_grid_view_(other.assembly_grid_view_)
    , source_space_(other.source_space_)
    , range_space_(other.range_space_)
    , symmetry_prefactor_(other.symmetry_prefactor_)
    , inner_penalty_(other.inner_penalty_)
    , dirichlet_penalty_(other.dirichlet_penalty_)
    , diffusion_(other.diffusion_)
    , weight_function_(other.weight_function_)
    , intersection_diameter_(other.intersection_diameter_)
    , element_mapper_(assembly_grid_view_)
  {}

  LaplaceIpdgFluxReconstructionOperator(ThisType&&) = default;

  bool linear() const override final
  {
    return true;
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  using BaseType::apply;

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param = {}) const override final
  {
    using D = typename AssemblyGridView::ctype;
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!source_space_.contains(source),
                  Exceptions::operator_error,
                  "source is not contained in source_space()!\n   source.size() = "
                      << source.size() << "\n   source_space().mappe().size() = " << source_space_.mapper().size());
    DUNE_THROW_IF(!range_space_.contains(range),
                  Exceptions::operator_error,
                  "range is not contained in range_space()!\n   range.size() = "
                      << range.size() << "\n   range_space().mappe().size() = " << range_space_.mapper().size());
    const auto source_function = make_discrete_function(source_space_, source);
    auto range_function = make_discrete_function(range_space_, range);
    // some preparations
    auto local_range = range_function.local_discrete_function();
    auto local_source_element = source_function.local_function();
    auto local_source_neighbor = source_function.local_function();
    auto local_diffusion_element = diffusion_.local_function();
    auto local_diffusion_neighbor = diffusion_.local_function();
    auto local_weight_element = weight_function_.local_function();
    auto local_weight_neighbor = weight_function_.local_function();
    auto rt_basis = range_space_.basis().localize();
    for (auto&& element : elements(assembly_grid_view_)) {
      local_range->bind(element);
      local_source_element->bind(element);
      local_diffusion_element->bind(element);
      local_weight_element->bind(element);
      rt_basis->bind(element);
      const auto& rt_fe = rt_basis->finite_element();
      // prepare
      const size_t sz = rt_basis->size();
      std::vector<bool> local_key_was_handled(sz, false);
      // determine the face dofs, therefore walk the intersections
      for (auto&& intersection : intersections(assembly_grid_view_, element)) {
        const auto intersection_to_local_key_map = rt_fe.coefficients().local_key_indices(1);
        const auto intersection_index = intersection.indexInInside();
        const auto& local_keys_assosiated_with_intersection = intersection_to_local_key_map[intersection_index];
        if (local_keys_assosiated_with_intersection.size() > 0) {
          const auto intersection_fe =
              make_local_orthonormal_finite_element<D, d - 1, F>(intersection.type(), rt_fe.order());
          const auto& intersection_Pk_basis = intersection_fe->basis();
          DUNE_THROW_IF(intersection_Pk_basis.size() != local_keys_assosiated_with_intersection.size(),
                        Exceptions::interpolation_error,
                        "intersection_Pk_basis.size() = " << intersection_Pk_basis.size()
                                                          << "\n   local_keys_assosiated_with_intersection.size() = "
                                                          << local_keys_assosiated_with_intersection.size());
          XT::LA::CommonDenseMatrix<F> lhs(
              local_keys_assosiated_with_intersection.size(), intersection_Pk_basis.size(), 0);
          XT::LA::CommonDenseVector<F> rhs(intersection_Pk_basis.size(), 0);
          bool there_are_intersection_dofs_to_determine = false;
          if (intersection.neighbor()) {
            const auto neighbor = intersection.outside();
            // only look at each intersection once
            if (element_mapper_.global_index(element, 0) < element_mapper_.global_index(neighbor, 0)) {
              there_are_intersection_dofs_to_determine = true;
              local_source_neighbor->bind(neighbor);
              local_diffusion_neighbor->bind(neighbor);
              local_weight_neighbor->bind(neighbor);
              // do a face quadrature
              const int max_polorder =
                  std::max(intersection_Pk_basis.order(),
                           std::max(local_source_element->order(param),
                                    std::max(intersection_Pk_basis.order(), local_source_neighbor->order(param))));
              for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                       intersection.type(), 2 * std::max(max_polorder, rt_basis->order()))) {
                const auto point_on_reference_intersection = quadrature_point.position();
                const auto point_in_reference_element =
                    intersection.geometryInInside().global(point_on_reference_intersection);
                const auto point_in_reference_neighbor =
                    intersection.geometryInOutside().global(point_on_reference_intersection);
                const auto quadrature_weight = quadrature_point.weight();
                const auto normal = intersection.unitOuterNormal(point_on_reference_intersection);
                const auto integration_factor =
                    intersection.geometry().integrationElement(point_on_reference_intersection);
                // evaluate bases ...
                const auto rt_basis_values = rt_basis->evaluate_set(point_in_reference_element);
                const auto intersection_Pk_basis_values =
                    intersection_Pk_basis.evaluate(point_on_reference_intersection);
                // ... and data functions
                const auto source_value_element = local_source_element->evaluate(point_in_reference_element, param);
                const auto source_value_neighbor = local_source_neighbor->evaluate(point_in_reference_neighbor, param);
                const auto source_grad_element = local_source_element->jacobian(point_in_reference_element, param)[0];
                const auto source_grad_neighbor =
                    local_source_neighbor->jacobian(point_in_reference_neighbor, param)[0];
                const auto diffusion_element = local_diffusion_element->evaluate(point_in_reference_element, param);
                const auto diffusion_neighbor = local_diffusion_neighbor->evaluate(point_in_reference_neighbor, param);
                // copied from LocalIPDGIntegrands::InnerPenalty and LocalLaplaceIPDGIntegrands::InnerCoupling
                const auto weight_element = local_weight_element->evaluate(point_in_reference_element, param);
                const auto weight_neighbor = local_weight_neighbor->evaluate(point_in_reference_neighbor, param);
                // compute the weighted penalty ...
                const double delta_plus = normal * (weight_neighbor * normal);
                const double delta_minus = normal * (weight_element * normal);
                const auto weight_minus = delta_plus / (delta_plus + delta_minus);
                const auto weight_plus = delta_minus / (delta_plus + delta_minus);
                const auto weight = (delta_plus * delta_minus) / (delta_plus + delta_minus); // half harmonic average
                const auto h = intersection_diameter_(intersection);
                const auto penalty = (inner_penalty_ * weight) / h;
                // ... and finally compute the LHS and RHS
                for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
                  const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
                  for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                    lhs.add_to_entry(ii,
                                     jj,
                                     quadrature_weight * integration_factor
                                         * (rt_basis_values[local_key_index] * normal)
                                         * intersection_Pk_basis_values[jj]);
                }
                for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                  rhs[jj] += quadrature_weight * integration_factor
                             * (penalty * (source_value_element - source_value_neighbor)
                                - symmetry_prefactor_
                                      * (normal
                                         * (weight_minus * (diffusion_element * source_grad_element)
                                            + weight_plus * (diffusion_neighbor * source_grad_neighbor))))
                             * intersection_Pk_basis_values[jj];
              }
            } else {
              there_are_intersection_dofs_to_determine = false;
              // we still need to mark the local DoFs as handled, since they were determined earlier when looking at
              // this intersection from the other side
              for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
                const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
                assert(!local_key_was_handled[local_key_index]);
                local_key_was_handled[local_key_index] = true;
              }
            }
          } else {
            there_are_intersection_dofs_to_determine = true;
            // do a face quadrature
            const int max_polorder = std::max(intersection_Pk_basis.order(), local_source_element->order(param));
            for (auto&& quadrature_point :
                 QuadratureRules<D, d - 1>::rule(intersection.type(), 2 * std::max(max_polorder, rt_basis->order()))) {
              const auto point_on_reference_intersection = quadrature_point.position();
              const auto point_in_reference_element =
                  intersection.geometryInInside().global(point_on_reference_intersection);
              const auto quadrature_weight = quadrature_point.weight();
              const auto normal = intersection.unitOuterNormal(point_on_reference_intersection);
              const auto integration_factor =
                  intersection.geometry().integrationElement(point_on_reference_intersection);
              // evaluate bases ...
              const auto rt_basis_values = rt_basis->evaluate_set(point_in_reference_element);
              const auto intersection_Pk_basis_values = intersection_Pk_basis.evaluate(point_on_reference_intersection);
              // ... and data functions ...
              const auto source_value_element = local_source_element->evaluate(point_in_reference_element, param);
              const auto source_grad_element = local_source_element->jacobian(point_in_reference_element, param)[0];
              const auto diffusion_element = local_diffusion_element->evaluate(point_in_reference_element, param);
              // copied from LocalIPDGIntegrands::BoundaryPenalty and LocalLaplaceIPDGIntegrands::DirichletCoupling
              const auto weight = local_weight_element->evaluate(point_in_reference_element, param);
              // compute the weighted penalty ...
              const auto h = intersection_diameter_(intersection);
              const auto penalty = (dirichlet_penalty_ * (normal * (weight * normal))) / h;
              // ... and finally compute the LHS and RHS
              for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
                const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
                for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                  lhs.add_to_entry(ii,
                                   jj,
                                   quadrature_weight * integration_factor * (rt_basis_values[local_key_index] * normal)
                                       * intersection_Pk_basis_values[jj]);
              }
              for (size_t jj = 0; jj < intersection_Pk_basis.size(); ++jj)
                rhs[jj] += quadrature_weight * integration_factor
                           * (penalty * source_value_element
                              - symmetry_prefactor_ * (normal * (diffusion_element * source_grad_element)))
                           * intersection_Pk_basis_values[jj];
            }
          }
          if (there_are_intersection_dofs_to_determine) {
            XT::LA::CommonDenseVector<F> intersection_dofs(local_keys_assosiated_with_intersection.size(), 0);
            try {
              intersection_dofs = XT::LA::solve(lhs, rhs);
            } catch (const XT::LA::Exceptions::linear_solver_failed& ee) {
              DUNE_THROW(Exceptions::interpolation_error,
                         "Failed to solve for DoFs associated with intersection "
                             << intersection_index << ", this was the original error:\n   " << ee.what());
            }
            for (size_t ii = 0; ii < local_keys_assosiated_with_intersection.size(); ++ii) {
              const size_t local_key_index = local_keys_assosiated_with_intersection[ii];
              assert(!local_key_was_handled[local_key_index]);
              local_range->dofs()[local_key_index] = intersection_dofs[ii];
              local_key_was_handled[local_key_index] = true;
            }
          }
        }
      }
      // determine the volume dofs
      // ... to be added ...
      // final checks that there are no other dofs left
      for (size_t ii = 0; ii < sz; ++ii)
        DUNE_THROW_IF(!local_key_was_handled[ii],
                      Exceptions::interpolation_error,
                      "The following DoF is neither associated with an intersection, nor with the element!"
                          << "\n   local DoF index: " << ii
                          << "\n   associated local_key: " << rt_fe.coefficients().local_key(ii));
    }
  } // ... apply(...)

private:
  const AssemblyGridViewType assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const double symmetry_prefactor_;
  const double inner_penalty_;
  const double dirichlet_penalty_;
  const XT::Functions::GridFunction<E, d, d> diffusion_;
  const XT::Functions::GridFunction<E, d, d> weight_function_;
  const std::function<double(const I&)> intersection_diameter_;
  const FiniteVolumeMapper<AssemblyGridViewType> element_mapper_;
}; // class LaplaceIpdgFluxReconstructionOperator


template <class M, class AGV, class SGV, class RGV>
LaplaceIpdgFluxReconstructionOperator<M, AGV, SGV, RGV> make_laplace_ipdg_flux_reconstruction_operator(
    AGV assembly_grid_view,
    const SpaceInterface<SGV>& source_space,
    const SpaceInterface<RGV, RGV::dimension>& range_space,
    const double& symmetry_prefactor,
    const double& inner_penalty,
    const double& dirichlet_penalty,
    XT::Functions::GridFunction<XT::Grid::extract_entity_t<AGV>, RGV::dimension, RGV::dimension> diffusion,
    XT::Functions::GridFunction<XT::Grid::extract_entity_t<AGV>, RGV::dimension, RGV::dimension> weight_function = {1.},
    const std::function<double(const XT::Grid::extract_intersection_t<AGV>&)>& intersection_diameter =
        LocalIPDGIntegrands::internal::default_intersection_diameter<XT::Grid::extract_intersection_t<AGV>>())
{
  return LaplaceIpdgFluxReconstructionOperator<M, AGV, SGV, RGV>(assembly_grid_view,
                                                                 source_space,
                                                                 range_space,
                                                                 symmetry_prefactor,
                                                                 inner_penalty,
                                                                 dirichlet_penalty,
                                                                 diffusion,
                                                                 weight_function,
                                                                 intersection_diameter);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LAPLACE_IPDG_FLUX_RECONSTRUCTION_HH
