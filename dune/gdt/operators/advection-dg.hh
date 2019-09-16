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

#ifndef DUNE_GDT_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_OPERATORS_ADVECTION_DG_HH

#include <cmath>

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/conversion.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/operators/advection-dg.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/local-mass-matrix.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


static inline double advection_dg_artificial_viscosity_default_nu_1()
{
  return 0.2;
}

static inline double advection_dg_artificial_viscosity_default_alpha_1()
{
  return 1.0;
}

static inline size_t advection_dg_artificial_viscosity_default_component()
{
  return 0;
}


/**
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class M, class AGV, size_t m = 1, class RGV = AGV, class SGV = AGV>
class AdvectionDgOperator : public LocalizableOperator<M, AGV, m, 1, m, 1, RGV, SGV>
{
  using ThisType = AdvectionDgOperator<M, AGV, m, RGV, SGV>;
  using BaseType = LocalizableOperator<M, AGV, m, 1, m, 1, RGV, SGV>;

protected:
  static const constexpr size_t d = SGV::dimension;
  using D = typename SGV::ctype;

public:
  using typename BaseType::F;
  using typename BaseType::V;

  using I = XT::Grid::extract_intersection_t<SGV>;
  using NumericalFluxType = NumericalFluxInterface<I, d, m, F>;
  using BoundaryTreatmentByCustomNumericalFluxOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomNumericalFluxOperator<I, V, SGV, m, F, F, RGV, V>;
  using BoundaryTreatmentByCustomExtrapolationOperatorType =
      LocalAdvectionDgBoundaryTreatmentByCustomExtrapolationOperator<I, V, SGV, m, F, F, RGV, V>;

  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  AdvectionDgOperator(
      const SGV& assembly_grid_view,
      const NumericalFluxType& numerical_flux,
      const SourceSpaceType& source_space,
      const RangeSpaceType& range_space,
      const XT::Grid::IntersectionFilter<SGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<SGV>(),
      const double& artificial_viscosity_nu_1 = advection_dg_artificial_viscosity_default_nu_1(),
      const double& artificial_viscosity_alpha_1 = advection_dg_artificial_viscosity_default_alpha_1(),
      const size_t artificial_viscosity_component = advection_dg_artificial_viscosity_default_component())
    : BaseType(assembly_grid_view, source_space, range_space)
    , numerical_flux_(numerical_flux.copy())
    , periodicity_exception_(periodicity_exception.copy())
    , local_mass_matrix_provider_(assembly_grid_view, range_space)
    , artificial_viscosity_nu_1_(artificial_viscosity_nu_1)
    , artificial_viscosity_alpha_1_(artificial_viscosity_alpha_1)
    , artificial_viscosity_component_(artificial_viscosity_component)
  {
    // we assemble these once, to be used in each apply later on
    auto walker = XT::Grid::make_walker(assembly_grid_view);
    walker.append(local_mass_matrix_provider_);
    walker.walk(/*use_tbb=*/true);
    // element contributions
    this->append(
        LocalAdvectionDgVolumeOperator<V, SGV, m, F, F, RGV, V>(local_mass_matrix_provider_, numerical_flux_->flux()));
    // contributions from inner intersections
    this->append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                     local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                 XT::Grid::ApplyOn::InnerIntersections<SGV>());
    // contributions from periodic boundaries
    this->append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                     local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                 *(XT::Grid::ApplyOn::PeriodicBoundaryIntersections<SGV>() && !(*periodicity_exception_)));
    // artificial viscosity by shock capturing [DF2015, Sec. 8.5]
    this->append(LocalAdvectionDgArtificialViscosityShockCapturingOperator<V, SGV, m, F, F, RGV, V>(
        local_mass_matrix_provider_,
        this->assembly_grid_view_,
        artificial_viscosity_nu_1_,
        artificial_viscosity_alpha_1_,
        artificial_viscosity_component_));
  } // AdvectionDgOperator(...)

  AdvectionDgOperator(ThisType&& source) = default;

  using BaseType::append;

  /// \name Non-periodic boundary treatment
  /// \{

  ThisType&
  append(typename BoundaryTreatmentByCustomNumericalFluxOperatorType::LambdaType numerical_boundary_treatment_flux,
         const int numerical_boundary_treatment_flux_order,
         const XT::Common::ParameterType& boundary_treatment_parameter_type = {},
         const XT::Grid::IntersectionFilter<SGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<SGV>())
  {
    this->append(BoundaryTreatmentByCustomNumericalFluxOperatorType(local_mass_matrix_provider_,
                                                                    numerical_boundary_treatment_flux,
                                                                    numerical_boundary_treatment_flux_order,
                                                                    boundary_treatment_parameter_type),
                 filter.copy());
    return *this;
  } // ... append(...)

  ThisType& append(typename BoundaryTreatmentByCustomExtrapolationOperatorType::LambdaType extrapolation,
                   const XT::Common::ParameterType& extrapolation_parameter_type = {},
                   const XT::Grid::IntersectionFilter<SGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<SGV>())
  {
    this->append(BoundaryTreatmentByCustomExtrapolationOperatorType(
                     local_mass_matrix_provider_, *numerical_flux_, extrapolation, extrapolation_parameter_type),
                 filter.copy());
    return *this;
  }

  /// \}

protected:
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  std::unique_ptr<XT::Grid::IntersectionFilter<SGV>> periodicity_exception_;
  LocalMassMatrixProvider<RGV, m, 1, F> local_mass_matrix_provider_;
  const double artificial_viscosity_nu_1_;
  const double artificial_viscosity_alpha_1_;
  const size_t artificial_viscosity_component_;
}; // class AdvectionDgOperator


template <class MatrixType, class AGV, size_t m, class F, class RGV, class SGV>
std::enable_if_t<XT::LA::is_matrix<MatrixType>::value, AdvectionDgOperator<MatrixType, AGV, m, RGV, SGV>>
make_advection_dg_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<XT::Grid::extract_intersection_t<AGV>, AGV::dimension, m, F>& numerical_flux,
    const SpaceInterface<SGV, m, 1, F>& source_space,
    const SpaceInterface<RGV, m, 1, F>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>(),
    const double& artificial_viscosity_nu_1 = advection_dg_artificial_viscosity_default_nu_1(),
    const double& artificial_viscosity_alpha_1 = advection_dg_artificial_viscosity_default_alpha_1(),
    const size_t artificial_viscosity_component = advection_dg_artificial_viscosity_default_component())
{
  return AdvectionDgOperator<MatrixType, AGV, m, RGV, SGV>(assembly_grid_view,
                                                           numerical_flux,
                                                           source_space,
                                                           range_space,
                                                           periodicity_exception,
                                                           artificial_viscosity_nu_1,
                                                           artificial_viscosity_alpha_1,
                                                           artificial_viscosity_component);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_DG_HH
