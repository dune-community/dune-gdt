// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_DG_HH
#define DUNE_GDT_OPERATORS_ADVECTION_DG_HH

#include <functional>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/local/operators/advection-dg.hh>
#include <dune/gdt/local/operators/advection-fv.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// forward
template <class DF, class GL>
class AdvectionDgOperator;


namespace internal {


template <class DF, class GL>
class AdvectionDgOperatorTraits
{
public:
  using derived_type = AdvectionDgOperator<DF, GL>;
  using FieldType = double;
  using JacobianType = int;
};


} // namespace internal


template <class DF, class GL = typename DF::SpaceType::GridLayerType>
class AdvectionDgOperator : public OperatorInterface<internal::AdvectionDgOperatorTraits<DF, GL>>
{
  static_assert(is_discrete_function<DF>::value, "");
  using SpaceType = typename DF::SpaceType;
  static_assert(XT::Grid::is_layer<GL>::value, "");
  using ThisType = AdvectionDgOperator<DF, GL>;
  using BaseType = OperatorInterface<internal::AdvectionDgOperatorTraits<DF, GL>>;
  using LocalVolumeOperatorType = LocalAdvectionDgVolumeOperator<SpaceType>;
  using LocalCouplingOperatorType = LocalAdvectionDgCouplingOperator<SpaceType>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::JacobianType;
  using NumericalFluxType = typename LocalCouplingOperatorType::NumericalFluxType;
  using LocalAdvectionDgBoundaryOperatorByCustomExtrapolationType =
      LocalAdvectionDgBoundaryOperatorByCustomExtrapolation<SpaceType>;
  using LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxType =
      LocalAdvectionDgBoundaryOperatorByCustomNumericalFlux<SpaceType>;

  AdvectionDgOperator(
      const GL& grid_layer,
      const NumericalFluxType& numerical_flx,
      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new XT::Grid::ApplyOn::NoIntersections<GL>())
    : BaseType(numerical_flx.parameter_type())
    , grid_layer_(grid_layer)
    , numerical_flux_(numerical_flx)
    , local_volume_operator_(numerical_flux_.access())
    , local_coupling_operator_(numerical_flux_.access())
    , periodicity_exception_(periodicity_exception)
  {
  }

  //  AdvectionFvOperator(
  //      const GL& grid_layer,
  //      NumericalFluxType*&& numerical_flx_ptr,
  //      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new
  //      XT::Grid::ApplyOn::NoIntersections<GL>())
  //    : BaseType(numerical_flx_ptr->parameter_type())
  //    , grid_layer_(grid_layer)
  //    , numerical_flux_(std::move(numerical_flx_ptr))
  //    , local_coupling_operator_(numerical_flux_.access())
  //    , periodicity_exception_(periodicity_exception)
  //  {
  //  }

  //  AdvectionFvOperator(
  //      const GL& grid_layer,
  //      const typename NumericalFluxType::FluxType& flux,
  //      std::function<typename NumericalFluxType::RangeType(const typename NumericalFluxType::RangeType&,
  //                                                          const typename NumericalFluxType::RangeType&,
  //                                                          const typename NumericalFluxType::DomainType&,
  //                                                          const XT::Common::Parameter&)> numerical_flux_lambda,
  //      const XT::Common::ParameterType& numerical_flux_parameter_type = {},
  //      XT::Grid::ApplyOn::WhichIntersection<GL>*&& periodicity_exception = new
  //      XT::Grid::ApplyOn::NoIntersections<GL>())
  //    : BaseType(numerical_flux_parameter_type)
  //    , grid_layer_(grid_layer)
  //    , numerical_flux_(new NumericalLambdaFlux<DF>(flux, numerical_flux_lambda, numerical_flux_parameter_type))
  //    , local_coupling_operator_(numerical_flux_.access())
  //    , periodicity_exception_(periodicity_exception)
  //  {
  //  }

  AdvectionDgOperator(const ThisType& other) = delete;
  AdvectionDgOperator(ThisType&& source) = delete;

  const NumericalFluxType& numerical_flux() const
  {
    return numerical_flux_.access();
  }

  void append(typename LocalAdvectionDgBoundaryOperatorByCustomExtrapolationType::LambdaType boundary_treatment_lambda,
              XT::Grid::ApplyOn::WhichIntersection<GL>*&& filter,
              const XT::Common::ParameterType& param_type = {})
  {
    boundary_treatment_by_extrapolation_operators_.emplace_back(
        new LocalAdvectionDgBoundaryOperatorByCustomExtrapolationType(
            numerical_flux(), boundary_treatment_lambda, param_type),
        std::move(filter));
  }

  void
  append(typename LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxType::LambdaType boundary_numerical_flux_lambda,
         const int flux_order,
         XT::Grid::ApplyOn::WhichIntersection<GL>*&& filter,
         const XT::Common::ParameterType& param_type = {})
  {
    boundary_treatment_by_boundary_flux_operators_.emplace_back(
        new LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxType(
            boundary_numerical_flux_lambda, flux_order, param_type),
        std::move(filter));
  }

  void apply(const DF& source, DF& range, const XT::Common::Parameter& param = {}) const
  {
    if (!source.vector().valid())
      DUNE_THROW(InvalidStateException, "source contains inf or nan!");
    range.vector() *= 0.;
    if (!range.vector().valid())
      DUNE_THROW(InvalidStateException, "range contains inf or nan!");
    LocalizableOperatorBase<GL, DF, DF> walker(grid_layer_, source, range);
    walker.append(local_volume_operator_, new XT::Grid::ApplyOn::AllEntities<GL>(), param, "dg_volume");
    walker.append(
        local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GL>(), param, "dg_coupling_inner");
    walker.append(local_coupling_operator_,
                  XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GL>() && !(*periodicity_exception_),
                  param,
                  "dg_coupling_periodic");
    for (const auto& boundary_treatment_op_and_filter : boundary_treatment_by_extrapolation_operators_) {
      const auto& op = boundary_treatment_op_and_filter.first.access();
      const auto& filter = *boundary_treatment_op_and_filter.second;
      walker.append(op, filter.copy(), param, "dg_boundary_extrapolation");
    }
    for (const auto& lambda_boundary_op_and_filter : boundary_treatment_by_boundary_flux_operators_) {
      const auto& op = lambda_boundary_op_and_filter.first.access();
      const auto& filter = *lambda_boundary_op_and_filter.second;
      walker.append(op, filter.copy(), param, "dg_boundary_custom_flux");
    }
    walker.walk();
    // apply the inverse mass matrix
    walker.append([&](const auto& entity) {
      using D = typename GL::ctype;
      static const size_t d = GL::dimension;
      auto local_range = range.local_discrete_function(entity);
      const auto& local_basis = local_range->basis();
      // prepare local containers
      XT::LA::CommonDenseMatrix<FieldType> local_mass_matrix(local_basis.size(), local_basis.size(), 0.);
      // create volume quadrature
      const auto volume_integrand_order = 2 * local_basis.order();
      // loop over all quadrature points
      for (const auto& quadrature_point :
           QuadratureRules<D, d>::rule(entity.type(), boost::numeric_cast<int>(volume_integrand_order))) {
        const auto xx = quadrature_point.position();
        // integration factors
        const auto integration_factor = entity.geometry().integrationElement(xx);
        const auto quadrature_weight = quadrature_point.weight();
        // evaluate everything
        const auto basis_values = local_basis.evaluate(xx);
        // compute mass matrix
        for (size_t ii = 0; ii < local_basis.size(); ++ii)
          for (size_t jj = 0; jj < local_basis.size(); ++jj)
            local_mass_matrix.add_to_entry(
                ii, jj, integration_factor * quadrature_weight * (basis_values[ii] * basis_values[jj]));
      }
      // solve
      XT::LA::CommonDenseVector<FieldType> tmp_rhs(local_basis.size(), 0.);
      XT::LA::CommonDenseVector<FieldType> tmp_solution(local_basis.size(), 0.);
      for (size_t ii = 0; ii < local_basis.size(); ++ii)
        tmp_rhs[ii] = local_range->vector().get(ii);
      XT::LA::make_solver(local_mass_matrix).apply(tmp_rhs, tmp_solution);
      for (size_t ii = 0; ii < local_basis.size(); ++ii)
        local_range->vector().set(ii, tmp_solution[ii]);
    });
    walker.walk();
    if (!range.vector().valid())
      DUNE_THROW(InvalidStateException, "range contains inf or nan!");
  } // ... apply(...)

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/,
                   const SourceType& /*source*/,
                   const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return FieldType();
  }

  template <class SourceType>
  JacobianType jacobian(const SourceType& /*source*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
    return JacobianType();
  }

  template <class SourceType>
  void
  jacobian(const SourceType& /*source*/, JacobianType& /*jac*/, const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/,
                     SourceType& /*source*/,
                     const XT::Common::Configuration& /*opts*/,
                     const Dune::XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(NotImplemented, "");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "");
    return std::vector<std::string>();
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "");
    return XT::Common::Configuration();
  }

private:
  const GL& grid_layer_;
  const XT::Common::ConstStorageProvider<NumericalFluxType> numerical_flux_;
  const LocalVolumeOperatorType local_volume_operator_;
  const LocalCouplingOperatorType local_coupling_operator_;
  std::list<std::pair<XT::Common::ConstStorageProvider<LocalAdvectionDgBoundaryOperatorByCustomExtrapolationType>,
                      std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>>>>
      boundary_treatment_by_extrapolation_operators_;
  std::list<std::pair<XT::Common::ConstStorageProvider<LocalAdvectionDgBoundaryOperatorByCustomNumericalFluxType>,
                      std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>>>>
      boundary_treatment_by_boundary_flux_operators_;
  std::unique_ptr<XT::Grid::ApplyOn::WhichIntersection<GL>> periodicity_exception_;
}; // class AdvectionDgOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_DG_HH
