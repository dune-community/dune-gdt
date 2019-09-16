// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_OPERATORS_ADVECTION_FV_HH

#include <dune/grid/common/partitionset.hh>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/local/assembler/operator-fd-jacobian-assemblers.hh>
#include <dune/gdt/local/operators/advection-fv.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


/**
 * \attention This operator will not work on a grid view with hanging nodes.
 *
 * \todo Refactor the coupling op as in the DG case to be applied on each side individually.
 *
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class M, class AGV, size_t m = 1, class RGV = AGV, class SGV = AGV>
class AdvectionFvOperator : public LocalizableOperator<M, AGV, m, 1, m, 1, RGV, SGV>
{
  using ThisType = AdvectionFvOperator<M, AGV, m, RGV, SGV>;
  using BaseType = LocalizableOperator<M, AGV, m, 1, m, 1, RGV, SGV>;

public:
  using BaseType::s_r;
  using BaseType::s_rC;
  using typename BaseType::F;
  using typename BaseType::V;

  using I = XT::Grid::extract_intersection_t<SGV>;
  using E = XT::Grid::extract_entity_t<SGV>;
  using NumericalFluxType = NumericalFluxInterface<I, SGV::dimension, m, F>;
  using BoundaryTreatmentByCustomNumericalFluxOperatorType =
      LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator<I, V, SGV, m, F, F, RGV, V>;
  using BoundaryTreatmentByCustomExtrapolationOperatorType =
      LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, V, SGV, m, F, F, RGV, V>;
  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, F>;

  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  AdvectionFvOperator(
      const SGV& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<SGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<SGV>())
    : BaseType(assembly_grid_view, source_space, range_space)
    , numerical_flux_(numerical_flux.copy())
    , periodicity_exception_(periodicity_exception.copy())
  {
    // contributions from inner intersections
    this->append(LocalAdvectionFvCouplingOperator<I, V, SGV, m, F, F, RGV, V>(*numerical_flux_),
                 XT::Grid::ApplyOn::InnerIntersectionsOnce<SGV>());
    // contributions from periodic boundaries
    this->append(LocalAdvectionFvCouplingOperator<I, V, SGV, m, F, F, RGV, V>(*numerical_flux_),
                 *(XT::Grid::ApplyOn::PeriodicBoundaryIntersectionsOnce<SGV>() && !(*periodicity_exception_)));
  }

  AdvectionFvOperator(ThisType&& source)
    : BaseType(std::move(source))
    , numerical_flux_(std::move(source.numerical_flux_))
    , periodicity_exception_(std::move(source.periodicity_exception_))
  {}

  using BaseType::append;

  /// \name Non-periodic boundary treatment
  /// \{

  ThisType&
  append(typename BoundaryTreatmentByCustomNumericalFluxOperatorType::LambdaType numerical_boundary_treatment_flux,
         const XT::Common::ParameterType& boundary_treatment_parameter_type = {},
         const XT::Grid::IntersectionFilter<SGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<SGV>())
  {
    this->append(BoundaryTreatmentByCustomNumericalFluxOperatorType(numerical_boundary_treatment_flux,
                                                                    boundary_treatment_parameter_type),
                 filter);
    return *this;
  }

  ThisType& append(typename BoundaryTreatmentByCustomExtrapolationOperatorType::LambdaType extrapolation,
                   const XT::Common::ParameterType& extrapolation_parameter_type = {},
                   const XT::Grid::IntersectionFilter<SGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<SGV>())
  {
    this->append(BoundaryTreatmentByCustomExtrapolationOperatorType(
                     *numerical_flux_, extrapolation, extrapolation_parameter_type),
                 filter);
    return *this;
  }

  /// \}

private:
  std::unique_ptr<const NumericalFluxType> numerical_flux_;
  std::unique_ptr<XT::Grid::IntersectionFilter<SGV>> periodicity_exception_;
}; // class AdvectionFvOperator


template <class MatrixType, class AGV, size_t m, class F, class RGV, class SGV>
std::enable_if_t<XT::LA::is_matrix<MatrixType>::value, AdvectionFvOperator<MatrixType, AGV, m, RGV, SGV>>
make_advection_fv_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<AGV>, AGV::dimension, m, F>& numerical_flux,
    const SpaceInterface<SGV, m, 1, F>& source_space,
    const SpaceInterface<RGV, m, 1, F>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>())
{
  return AdvectionFvOperator<MatrixType, AGV, m, RGV, SGV>(
      assembly_grid_view, numerical_flux, source_space, range_space, periodicity_exception);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_FV_HH
