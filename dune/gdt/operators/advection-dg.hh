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
template <class M, class SGV, size_t m = 1, class RGV = SGV>
class AdvectionDgOperator : public OperatorInterface<M, SGV, m, 1, m, 1, RGV>
{
  using ThisType = AdvectionDgOperator<M, SGV, m, RGV>;
  using BaseType = OperatorInterface<M, SGV, m, 1, m, 1, RGV>;

protected:
  static const constexpr size_t d = SGV::dimension;
  using D = typename SGV::ctype;

public:
  using typename BaseType::F;
  using typename BaseType::V;

  using NumericalFluxType = NumericalFluxInterface<d, m, F>;

  using I = XT::Grid::extract_intersection_t<SGV>;
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
    : BaseType(numerical_flux.parameter_type())
    , assembly_grid_view_(assembly_grid_view)
    , numerical_flux_(numerical_flux.copy())
    , source_space_(source_space)
    , range_space_(range_space)
    , periodicity_exception_(periodicity_exception.copy())
    , local_mass_matrix_provider_(assembly_grid_view_, range_space_)
    , artificial_viscosity_nu_1_(artificial_viscosity_nu_1)
    , artificial_viscosity_alpha_1_(artificial_viscosity_alpha_1)
    , artificial_viscosity_component_(artificial_viscosity_component)
  {
    // we assemble these once, to be used in each apply later on
    auto walker = XT::Grid::make_walker(assembly_grid_view_);
    walker.append(local_mass_matrix_provider_);
    walker.walk(/*use_tbb=*/true);
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
         const XT::Grid::IntersectionFilter<SGV>& filter = XT::Grid::ApplyOn::BoundaryIntersections<SGV>())
  {
    boundary_treatments_by_custom_numerical_flux_.emplace_back(
        new BoundaryTreatmentByCustomNumericalFluxOperatorType(local_mass_matrix_provider_,
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
    boundary_treatments_by_custom_extrapolation_.emplace_back(
        new BoundaryTreatmentByCustomExtrapolationOperatorType(
            local_mass_matrix_provider_, *numerical_flux_, extrapolation, extrapolation_parameter_type),
        filter.copy());
    return *this;
  }

  /// \}

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param = {}) const override
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
    auto source_function = make_discrete_function(this->source_space_, source);
    auto range_function = make_discrete_function(this->range_space_, range);
    // set up the actual operator
    auto localizable_op = make_localizable_operator(this->assembly_grid_view_, source_function, range_function);
    // element contributions
    localizable_op.append(
        LocalAdvectionDgVolumeOperator<V, SGV, m, F, F, RGV, V>(local_mass_matrix_provider_, numerical_flux_->flux()),
        param);
    // contributions from inner intersections
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                              local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                          param,
                          XT::Grid::ApplyOn::InnerIntersections<SGV>());
    // contributions from periodic boundaries
    localizable_op.append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                              local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                          param,
                          *(XT::Grid::ApplyOn::PeriodicBoundaryIntersections<SGV>() && !(*periodicity_exception_)));
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
    // artificial viscosity by shock capturing [DF2015, Sec. 8.5]
    localizable_op.append(LocalAdvectionDgArtificialViscosityShockCapturingOperator<V, SGV, m, F, F, RGV, V>(
                              local_mass_matrix_provider_,
                              this->assembly_grid_view_,
                              artificial_viscosity_nu_1_,
                              artificial_viscosity_alpha_1_,
                              artificial_viscosity_component_),
                          param);
    // and apply it in a grid walk
    localizable_op.assemble(/*use_tbb=*/true);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

  std::vector<std::string> jacobian_options() const override final
  {
    return {"finite-differences"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->jacobian_options().at(0), Exceptions::operator_error, "type = " << type);
    return {{"type", type}, {"eps", "1e-7"}};
  }

  using BaseType::jacobian;

  void jacobian(const VectorType& source,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, opts);
    DUNE_THROW_IF(opts.get<std::string>("type") != jacobian_options().at(0), Exceptions::operator_error, opts);
    const auto default_opts = jacobian_options(jacobian_options().at(0));
    const auto eps = opts.get("eps", default_opts.template get<double>("eps"));
    const auto parameter = param + XT::Common::Parameter({"finite-difference-jacobians.eps", eps});
    // append the same local ops with the same filters as in apply() above
    // element contributions
    jacobian_op.append(
        LocalAdvectionDgVolumeOperator<V, SGV, m, F, F, RGV, V>(local_mass_matrix_provider_, numerical_flux_->flux()),
        source,
        parameter);
    // contributions from inner intersections
    jacobian_op.append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                           local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                       source,
                       parameter,
                       XT::Grid::ApplyOn::InnerIntersections<SGV>());
    // contributions from periodic boundaries
    jacobian_op.append(LocalAdvectionDgCouplingOperator<I, V, SGV, m, F, F, RGV, V>(
                           local_mass_matrix_provider_, *numerical_flux_, /*compute_outside=*/false),
                       source,
                       parameter,
                       *(XT::Grid::ApplyOn::PeriodicBoundaryIntersections<SGV>() && !(*periodicity_exception_)));
    // contributions from other boundaries by custom numerical flux
    for (const auto& boundary_treatment : boundary_treatments_by_custom_numerical_flux_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      jacobian_op.append(boundary_op, source, parameter, filter);
    }
    // contributions from other boundaries by custom extrapolation
    for (const auto& boundary_treatment : boundary_treatments_by_custom_extrapolation_) {
      const auto& boundary_op = *boundary_treatment.first;
      const auto& filter = *boundary_treatment.second;
      jacobian_op.append(boundary_op, source, parameter, filter);
    }
    // artificial viscosity by shock capturing [DF2015, Sec. 8.5]
    jacobian_op.append(LocalAdvectionDgArtificialViscosityShockCapturingOperator<V, SGV, m, F, F, RGV, V>(
                           local_mass_matrix_provider_,
                           this->assembly_grid_view_,
                           artificial_viscosity_nu_1_,
                           artificial_viscosity_alpha_1_,
                           artificial_viscosity_component_),
                       source,
                       param);
  } // ... jacobian(...)

protected:
  const SGV assembly_grid_view_;
  const std::unique_ptr<const NumericalFluxType> numerical_flux_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::unique_ptr<XT::Grid::IntersectionFilter<SGV>> periodicity_exception_;
  LocalMassMatrixProvider<RGV, m, 1, F> local_mass_matrix_provider_;
  const double artificial_viscosity_nu_1_;
  const double artificial_viscosity_alpha_1_;
  const size_t artificial_viscosity_component_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomNumericalFluxOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<SGV>>>>
      boundary_treatments_by_custom_numerical_flux_;
  std::list<std::pair<std::unique_ptr<BoundaryTreatmentByCustomExtrapolationOperatorType>,
                      std::unique_ptr<XT::Grid::IntersectionFilter<SGV>>>>
      boundary_treatments_by_custom_extrapolation_;
}; // class AdvectionDgOperator


template <class MatrixType, class SGV, size_t m, class F, class RGV>
std::enable_if_t<XT::LA::is_matrix<MatrixType>::value, AdvectionDgOperator<MatrixType, SGV, m, RGV>>
make_advection_dg_operator(
    const SGV& assembly_grid_view,
    const NumericalFluxInterface<SGV::dimension, m, F>& numerical_flux,
    const SpaceInterface<SGV, m, 1, F>& source_space,
    const SpaceInterface<RGV, m, 1, F>& range_space,
    const XT::Grid::IntersectionFilter<SGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<SGV>(),
    const double& artificial_viscosity_nu_1 = advection_dg_artificial_viscosity_default_nu_1(),
    const double& artificial_viscosity_alpha_1 = advection_dg_artificial_viscosity_default_alpha_1(),
    const size_t artificial_viscosity_component = advection_dg_artificial_viscosity_default_component())
{
  return AdvectionDgOperator<MatrixType, SGV, m, RGV>(assembly_grid_view,
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
