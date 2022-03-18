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
#include "operator.hh"

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
template <class AGV,
          size_t m = 1,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class RGV = AGV,
          class SGV = AGV>
class AdvectionFvOperator : public Operator<AGV, m, 1, m, 1, F, M, RGV, SGV>
{
public:
  using ThisType = AdvectionFvOperator;
  using BaseType = Operator<AGV, m, 1, m, 1, F, M, RGV, SGV>;

  using BaseType::s_r;
  using BaseType::s_rC;
  using typename BaseType::V;

  using I = XT::Grid::extract_intersection_t<AGV>;
  using E = XT::Grid::extract_entity_t<AGV>;
  using NumericalFluxType = NumericalFluxInterface<I, AGV::dimension, m, F>;
  using BoundaryTreatmentByCustomNumericalFluxOperatorType =
      LocalAdvectionFvBoundaryTreatmentByCustomNumericalFluxOperator<I, V, AGV, m, F, F, RGV, V>;
  using BoundaryTreatmentByCustomExtrapolationOperatorType =
      LocalAdvectionFvBoundaryTreatmentByCustomExtrapolationOperator<I, V, AGV, m, F, F, RGV, V>;

  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  AdvectionFvOperator(
      const AGV& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
      const std::string& logging_prefix = "",
      const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType(assembly_grid_view,
               source_space,
               range_space,
               /*requires_assembly_=*/false,
               logging_prefix.empty() ? "AdvectionFvOperator" : logging_prefix,
               logging_state)
    , numerical_flux_(numerical_flux.copy())
    , periodicity_exception_(periodicity_exception.copy())
  {
    // contributions from inner intersections
    *this += {LocalAdvectionFvCouplingOperator<I, V, AGV, m, F, F, RGV, V>(*numerical_flux_),
              XT::Grid::ApplyOn::InnerIntersectionsOnce<AGV>()};
    // contributions from periodic boundaries
    *this += {LocalAdvectionFvCouplingOperator<I, V, AGV, m, F, F, RGV, V>(*numerical_flux_),
              *(XT::Grid::ApplyOn::PeriodicBoundaryIntersectionsOnce<AGV>() && !(*periodicity_exception_))};
  }

  AdvectionFvOperator(ThisType&& source)
    : BaseType(std::move(source))
    , numerical_flux_(std::move(source.numerical_flux_))
    , periodicity_exception_(std::move(source.periodicity_exception_))
  {}

  /// \name These methods can be used to define non-periodic boundary treatment
  /// \{

  ThisType& boundary_treatment(
      typename BoundaryTreatmentByCustomNumericalFluxOperatorType::LambdaType numerical_boundary_treatment_flux,
      const XT::Common::ParameterType& boundary_treatment_parameter_type = {},
      const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    *this += {BoundaryTreatmentByCustomNumericalFluxOperatorType(numerical_boundary_treatment_flux,
                                                                 boundary_treatment_parameter_type),
              filter};
    return *this;
  }

  ThisType&
  boundary_treatment(typename BoundaryTreatmentByCustomExtrapolationOperatorType::LambdaType extrapolation,
                     const XT::Common::ParameterType& extrapolation_parameter_type = {},
                     const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<AGV>())
  {
    *this += {BoundaryTreatmentByCustomExtrapolationOperatorType(
                  *numerical_flux_, extrapolation, extrapolation_parameter_type),
              filter};
    return *this;
  }

  const NumericalFluxType& numerical_flux() const
  {
    return *numerical_flux_;
  }

  /// \}

private:
  std::unique_ptr<const NumericalFluxType> numerical_flux_;
  std::unique_ptr<XT::Grid::IntersectionFilter<AGV>> periodicity_exception_;
}; // class AdvectionFvOperator


template <class MatrixType, // <- has to be specified manually
          class AGV,
          size_t m,
          class F,
          class RGV,
          class SGV>
auto make_advection_fv_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<AGV>, AGV::dimension, m, F>& numerical_flux,
    const SpaceInterface<SGV, m, 1, F>& source_space,
    const SpaceInterface<RGV, m, 1, F>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
    const std::string& logging_prefix = "",
    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return AdvectionFvOperator<AGV, m, F, MatrixType, RGV, SGV>(assembly_grid_view,
                                                              numerical_flux,
                                                              source_space,
                                                              range_space,
                                                              periodicity_exception,
                                                              logging_prefix,
                                                              logging_state);
}

template <class AGV, size_t m, class F, class RGV, class SGV>
auto make_advection_fv_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<AGV>, AGV::dimension, m, F>& numerical_flux,
    const SpaceInterface<SGV, m, 1, F>& source_space,
    const SpaceInterface<RGV, m, 1, F>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
    const std::string& logging_prefix = "",
    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_advection_fv_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(assembly_grid_view,
                                                                         numerical_flux,
                                                                         source_space,
                                                                         range_space,
                                                                         periodicity_exception,
                                                                         logging_prefix,
                                                                         logging_state);
}


template <class MatrixType, // <- has to be specified manually
          class GV,
          size_t m,
          class F>
auto make_advection_fv_operator(
    const SpaceInterface<GV, m, 1, F>& space,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<GV>, GV::dimension, m, F>& numerical_flux,
    const XT::Grid::IntersectionFilter<GV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<GV>(),
    const std::string& logging_prefix = "",
    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return AdvectionFvOperator<GV, m, F, MatrixType, GV, GV>(
      space.grid_view(), numerical_flux, space, space, periodicity_exception, logging_prefix, logging_state);
}

template <class GV, size_t m, class F>
auto make_advection_fv_operator(
    const SpaceInterface<GV, m, 1, F>& space,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<GV>, GV::dimension, m, F>& numerical_flux,
    const XT::Grid::IntersectionFilter<GV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<GV>(),
    const std::string& logging_prefix = "",
    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_advection_fv_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(
      space, numerical_flux, periodicity_exception, logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_FV_HH
