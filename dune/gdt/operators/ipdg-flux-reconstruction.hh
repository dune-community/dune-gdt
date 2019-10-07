// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_OPERATORS_IPDG_FLUX_RECONSTRUCTION_HH
#define DUNE_GDT_OPERATORS_IPDG_FLUX_RECONSTRUCTION_HH

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/orthonormal.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/mapper/finite-volume.hh>

namespace Dune {
namespace GDT {


/**
 * \todo This is just the old implementation copied over. Update analoggously to the RT interpolation!
 */
template <class M,
          class AssemblyGridView,
          LocalEllipticIpdgIntegrands::Method ipdg = LocalEllipticIpdgIntegrands::default_method,
          class SGV = AssemblyGridView,
          class RGV = AssemblyGridView>
class IpdgFluxReconstructionOperator : public OperatorInterface<M, SGV, 1, 1, RGV::dimension, 1, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  using BaseType = OperatorInterface<M, SGV, 1, 1, RGV::dimension, 1, RGV>;
  using ThisType = IpdgFluxReconstructionOperator;

public:
  using typename BaseType::F;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AGV = AssemblyGridView;
  using AssemblyGridViewType = AGV;
  using E = XT::Grid::extract_entity_t<AssemblyGridViewType>;
  static const constexpr size_t d = RGV::dimension;

  IpdgFluxReconstructionOperator(AssemblyGridViewType assembly_grid_view,
                                 const SourceSpaceType& src_spc,
                                 const RangeSpaceType& rng_spc,
                                 const XT::Functions::GridFunctionInterface<E>& diffusion_factor,
                                 const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor)
    : assembly_grid_view_(assembly_grid_view)
    , source_space_(src_spc)
    , range_space_(rng_spc)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , element_mapper_(assembly_grid_view_)
  {
    DUNE_THROW_IF(range_space_.type() != SpaceType::raviart_thomas, Exceptions::operator_error, "");
    DUNE_THROW_IF(range_space_.max_polorder() != 0, Exceptions::operator_error, "Not implemented yet!");
  }

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
    using I = XT::Grid::extract_intersection_t<AssemblyGridViewType>;
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
    auto local_df_element = diffusion_factor_.local_function();
    auto local_df_neighbor = diffusion_factor_.local_function();
    auto local_dt_element = diffusion_tensor_.local_function();
    auto local_dt_neighbor = diffusion_tensor_.local_function();
    auto rt_basis = range_space_.basis().localize();
    for (auto&& element : elements(assembly_grid_view_)) {
      local_range->bind(element);
      local_source_element->bind(element);
      local_df_element->bind(element);
      local_dt_element->bind(element);
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
              make_local_orthonormal_finite_element<D, d - 1, F>(intersection.geometry().type(), rt_fe.order());
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
              local_df_neighbor->bind(neighbor);
              local_dt_neighbor->bind(neighbor);
              // do a face quadrature
              const int max_polorder =
                  std::max(intersection_Pk_basis.order(),
                           std::max(local_source_element->order(param),
                                    std::max(intersection_Pk_basis.order(), local_source_neighbor->order(param))));
              for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                       intersection.geometry().type(), 2 * std::max(max_polorder, rt_basis->order()))) {
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
                const auto df_value_element = local_df_element->evaluate(point_in_reference_element, param);
                const auto df_value_neighbor = local_df_neighbor->evaluate(point_in_reference_neighbor, param);
                const auto dt_value_element = local_dt_element->evaluate(point_in_reference_element, param);
                const auto dt_value_neighbor = local_dt_neighbor->evaluate(point_in_reference_neighbor, param);
                const auto diffusion_element = dt_value_element * df_value_element;
                const auto diffusion_neighbor = dt_value_neighbor * df_value_neighbor;
                // compute penalty factor (see Epshteyn, Riviere, 2007)
                const F sigma = LocalEllipticIpdgIntegrands::internal::inner_sigma(max_polorder);
                const double beta = LocalEllipticIpdgIntegrands::internal::default_beta(d);
                // compute weighting (see Ern, Stephansen, Zunino 2007)
                using IpdgHelper = typename LocalEllipticIpdgIntegrands::Inner<I, F, ipdg>::template IPDG<ipdg>;
                const F delta_plus =
                    IpdgHelper::delta_plus(df_value_neighbor, dt_value_neighbor, diffusion_neighbor, normal);
                const F delta_minus =
                    IpdgHelper::delta_minus(df_value_element, dt_value_element, diffusion_element, normal);
                const F gamma = IpdgHelper::gamma(delta_plus, delta_minus);
                const F penalty = IpdgHelper::penalty(df_value_element,
                                                      dt_value_neighbor,
                                                      df_value_neighbor,
                                                      dt_value_element,
                                                      normal,
                                                      sigma,
                                                      gamma,
                                                      intersection.geometry().volume(),
                                                      beta);
                const F weight_plus = IpdgHelper::weight_plus(delta_plus, delta_minus);
                const F weight_minus = IpdgHelper::weight_minus(delta_plus, delta_minus);
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
                                - normal
                                      * (weight_minus * (diffusion_element * source_grad_element)
                                         + weight_plus * (diffusion_neighbor * source_grad_neighbor)))
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
            for (auto&& quadrature_point : QuadratureRules<D, d - 1>::rule(
                     intersection.geometry().type(), 2 * std::max(max_polorder, rt_basis->order()))) {
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
              const auto df_value_element = local_df_element->evaluate(point_in_reference_element, param);
              const auto dt_value_element = local_dt_element->evaluate(point_in_reference_element, param);
              const auto diffusion_element = dt_value_element * df_value_element;
              // compute penalty (see Epshteyn, Riviere, 2007)
              const F sigma = LocalEllipticIpdgIntegrands::internal::boundary_sigma(max_polorder);
              const double beta = LocalEllipticIpdgIntegrands::internal::default_beta(d);
              // compute weighting (see Ern, Stephansen, Zunino 2007)
              using IpdgHelper =
                  typename LocalEllipticIpdgIntegrands::DirichletBoundaryLhs<I, F, ipdg>::template IPDG<ipdg>;
              const F gamma = IpdgHelper::gamma(diffusion_element, normal);
              const F penalty = IpdgHelper::penalty(sigma, gamma, intersection.geometry().volume(), beta);
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
                           * (penalty * source_value_element - normal * (diffusion_element * source_grad_element))
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
  const XT::Functions::GridFunctionInterface<E>& diffusion_factor_;
  const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor_;
  const FiniteVolumeMapper<AssemblyGridViewType> element_mapper_;
}; // class IpdgFluxReconstructionOperator


template <class MatrixType,
          class AssemblyGridViewType,
          class SGV,
          class RGV,
          LocalEllipticIpdgIntegrands::Method ipdg = LocalEllipticIpdgIntegrands::default_method>
IpdgFluxReconstructionOperator<MatrixType, AssemblyGridViewType, ipdg, SGV, RGV> make_ipdg_flux_reconstruction_operator(
    const AssemblyGridViewType& assembly_grid_view,
    const SpaceInterface<SGV>& source_space,
    const SpaceInterface<RGV, RGV::dimension>& range_space,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<AssemblyGridViewType>>& diffusion_factor,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<AssemblyGridViewType>,
                                               RGV::dimension,
                                               RGV::dimension>& diffusion_tensor)
{
  return IpdgFluxReconstructionOperator<MatrixType, AssemblyGridViewType, ipdg, SGV, RGV>(
      assembly_grid_view, source_space, range_space, diffusion_factor, diffusion_tensor);
}


template <class MatrixType, LocalEllipticIpdgIntegrands::Method ipdg, class AssemblyGridViewType, class SGV, class RGV>
IpdgFluxReconstructionOperator<MatrixType, AssemblyGridViewType, ipdg, SGV, RGV> make_ipdg_flux_reconstruction_operator(
    const AssemblyGridViewType& assembly_grid_view,
    const SpaceInterface<SGV>& source_space,
    const SpaceInterface<RGV, RGV::dimension>& range_space,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<AssemblyGridViewType>>& diffusion_factor,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<AssemblyGridViewType>,
                                               RGV::dimension,
                                               RGV::dimension>& diffusion_tensor)
{
  return IpdgFluxReconstructionOperator<MatrixType, AssemblyGridViewType, ipdg, SGV, RGV>(
      assembly_grid_view, source_space, range_space, diffusion_factor, diffusion_tensor);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_IPDG_FLUX_RECONSTRUCTION_HH
