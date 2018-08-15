// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_OPERATORS_ADVECTION_DG_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker/filters.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/advection-dg.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


/**
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class AssemblyGridView,
          class SV,
          size_t m = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          class RF = SF,
          class RGV = AssemblyGridView,
          class RV = SV>
class AdvectionDgOperator : public OperatorInterface<SV,
                                                     SGV,
                                                     m,
                                                     1,
                                                     SF,
                                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                                     m,
                                                     1,
                                                     RF,
                                                     RGV,
                                                     RV>
{
  // No need to check the rest, is done in OperatorInterface.
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert((AssemblyGridView::dimension == SGV::dimension) && (AssemblyGridView::dimension == RGV::dimension), "");

  using ThisType = AdvectionDgOperator<AssemblyGridView, SV, m, SF, SGV, RF, RGV, RV>;
  using BaseType = OperatorInterface<SV,
                                     SGV,
                                     m,
                                     1,
                                     SF,
                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                     m,
                                     1,
                                     RF,
                                     RGV,
                                     RV>;
  using AGV = AssemblyGridView;

public:
  using AssemblyGridViewType = AssemblyGridView;
  using NumericalFluxType = NumericalFluxInterface<AssemblyGridView::dimension, m, RF>;

  using I = XT::Grid::extract_intersection_t<AGV>;
  using BoundaryTreatmentByCustomNumericalFluxOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator<I, SV, SGV, m, SF, RF, RGV, RV>;
  using BoundaryTreatmentByCustomExtrapolationOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator<I, SV, SGV, m, SF, RF, RGV, RV>;

  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;

  using typename BaseType::SourceVectorType;
  using typename BaseType::RangeVectorType;

  AdvectionDgOperator(
      const AssemblyGridViewType& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>())
    : BaseType(numerical_flux.parameter_type())
    , assembly_grid_view_(assembly_grid_view)
    , numerical_flux_(numerical_flux.copy())
    , source_space_(source_space)
    , range_space_(range_space)
    , periodicity_exception_(periodicity_exception.copy())
  {
  }

  AdvectionDgOperator(ThisType&& source) = default;

  bool linear() const override final
  {
    return numerical_flux_->linear();
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  /// \name Non-periodic boundary treatment
  /// \{

  ThisType&
  append(typename BoundaryTreatmentByCustomNumericalFluxOperatorType::LambdaType numerical_boundary_treatment_flux,
         const int numerical_boundary_treatment_flux_order,
         const XT::Common::ParameterType& boundary_treatment_parameter_type = {},
         const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    boundary_treatments_by_custom_numerical_flux_.emplace_back(
        new BoundaryTreatmentByCustomNumericalFluxOperatorType(numerical_boundary_treatment_flux,
                                                               numerical_boundary_treatment_flux_order,
                                                               boundary_treatment_parameter_type),
        filter.copy());
    return *this;
  } // ... append(...)

  ThisType& append(typename BoundaryTreatmentByCustomExtrapolationOperatorType::LambdaType extrapolation,
                   const XT::Common::ParameterType& extrapolation_parameter_type = {},
                   const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    boundary_treatments_by_custom_extrapolation_.emplace_back(
        new BoundaryTreatmentByCustomExtrapolationOperatorType(
            *numerical_flux_, extrapolation, extrapolation_parameter_type),
        filter.copy());
    return *this;
  }

  /// \}

  using BaseType::apply;

  void apply(const SourceVectorType& source,
             RangeVectorType& range,
             const XT::Common::Parameter& param = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
    const auto source_function = make_discrete_function(source_space_, source);
    auto range_function = make_discrete_function(range_space_, range);
    // set up the actual operator
    auto localizable_op = make_localizable_operator(assembly_grid_view_, source_function, range_function);
    // element contributions
    localizable_op.append(LocalAdvectionDgVolumeOperator<SV, SGV, m, SF, RF, RGV, RV>(numerical_flux_->flux()), param);
    // contributions from inner intersections
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(*numerical_flux_),
                          param,
                          XT::Grid::ApplyOn::InnerIntersectionsOnce<AGV>());
    // contributions from periodic boundaries
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(*numerical_flux_),
                          param,
                          *(XT::Grid::ApplyOn::PeriodicBoundaryIntersectionsOnce<AGV>() && !(*periodicity_exception_)));
    // contributions from other boundaries by custom numerical flux
    for (const auto& boundary_treatment : boundary_treatments_by_custom_numerical_flux_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      localizable_op.append(boundary_op, param, filter);
    }
    // contributions from other boundaries by custom extrapolation
    for (const auto& boundary_treatment : boundary_treatments_by_custom_extrapolation_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      localizable_op.append(boundary_op, param, filter);
    }
    // do the actual (first) grid walk: the above operators will be applied and afterwards cleared
    localizable_op.assemble(/*use_tbb=*/true);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
    // apply the inverse mass matrix by appending a grid element functor, use the localizable_op as a XT::Grid::Walker
    localizable_op.append([&](const auto& element) {
      // prepare
      // (creating these objects before the grid walk and reusing them would be more efficient, but not thread safe)
      using E = XT::Grid::extract_entity_t<RGV>;
      const LocalElementIntegralBilinearForm<E, m, 1, RF, RF> local_l2_bilinear_form(
          LocalElementProductIntegrand<E, m, 1, RF, RF>(1.));
      auto local_range = range_function.local_discrete_function(element);
      auto local_rhs = XT::LA::convert_to<XT::LA::CommonDenseVector<RF>>(local_range->dofs());
      XT::LA::CommonDenseVector<RF> local_solution(range_space_.mapper().max_local_size(), 0.);
      // compute local L_2 matrix
      const auto& range_basis = local_range->basis();
      auto local_mass_matrix =
          XT::LA::convert_to<XT::LA::CommonDenseMatrix<RF>>(local_l2_bilinear_form.apply2(range_basis, range_basis));
      // solve
      XT::LA::solve(local_mass_matrix, local_rhs, local_solution);
      // set DoFs
      for (size_t ii = 0; ii < range_basis.size(); ++ii)
        local_range->dofs()[ii] = local_solution[ii];
    });
    // do the actual (second) grid walk
    localizable_op.walk(/*use_tbb=*/true); // Do not call assemble(), will do nothing the 2nd time!
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

private:
  const AssemblyGridViewType& assembly_grid_view_;
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::unique_ptr<XT::Grid::IntersectionFilter<AssemblyGridViewType>> periodicity_exception_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomNumericalFluxOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<AGV>>>>
      boundary_treatments_by_custom_numerical_flux_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomExtrapolationOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<AGV>>>>
      boundary_treatments_by_custom_extrapolation_;
}; // class AdvectionDgOperator


template <class VectorType, class AGV, size_t m, class RF, class SGV, class SF, class RGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value,
                 AdvectionDgOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>>
make_advection_dg_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<AGV::dimension, m, RF>& numerical_flux,
    const SpaceInterface<SGV, m, 1, SF>& source_space,
    const SpaceInterface<RGV, m, 1, RF>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>())
{
  return AdvectionDgOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>(
      assembly_grid_view, numerical_flux, source_space, range_space, periodicity_exception);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_DG_HH
