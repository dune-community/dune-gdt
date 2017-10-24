// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_OPERATORS_ADVECTION_FV_HH

#include <functional>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/type_traits.hh>
#include <dune/gdt/local/operators/advection-fv.hh>

#include "base.hh"

namespace Dune {
namespace GDT {


// forward
template <class DF, class GL>
class AdvectionFvOperator;


namespace internal {


template <class DF, class GL>
class AdvectionFvOperatorTraits
{
public:
  using derived_type = AdvectionFvOperator<DF, GL>;
  using FieldType = double;
  using JacobianType = int;
};


} // namespace internal


template <class DF, class GL = typename DF::SpaceType::GridLayerType>
class AdvectionFvOperator : public OperatorInterface<internal::AdvectionFvOperatorTraits<DF, GL>>
{
  static_assert(is_discrete_function<DF>::value, "");
  using SpaceType = typename DF::SpaceType;
  static_assert(SpaceType::polOrder == 0, "Use AdvectionDgOperator instead!");
  static_assert(XT::Grid::is_layer<GL>::value, "");
  using ThisType = AdvectionFvOperator<DF, GL>;
  using BaseType = OperatorInterface<internal::AdvectionFvOperatorTraits<DF, GL>>;
  using LocalCouplingOperatorType = LocalAdvectionFvInnerOperator<SpaceType>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::JacobianType;
  using NumericalFluxType = typename LocalCouplingOperatorType::NumericalFluxType;

  AdvectionFvOperator(const GL& grid_layer, const NumericalFluxType& numerical_flx)
    : grid_layer_(grid_layer)
    , numerical_flux_()
    , local_coupling_operator_(numerical_flux_.access())
  {
  }

  AdvectionFvOperator(
      const GL& grid_layer,
      const typename NumericalFluxType::FluxType& flux,
      std::function<typename NumericalFluxType::RangeType(const typename NumericalFluxType::RangeType&,
                                                          const typename NumericalFluxType::RangeType&,
                                                          const typename NumericalFluxType::DomainType&,
                                                          const XT::Common::Parameter&)> numerical_flux_lambda,
      const XT::Common::ParameterType& numerical_flux_parameter_type = {})
    : grid_layer_(grid_layer)
    , numerical_flux_(new NumericalLambdaFlux<DF>(flux, numerical_flux_lambda, numerical_flux_parameter_type))
    , local_coupling_operator_(numerical_flux_.access())
  {
  }

  AdvectionFvOperator(const ThisType& other) = delete;
  AdvectionFvOperator(ThisType&& source) = delete;

  const NumericalFluxType& numerical_flux() const
  {
    return numerical_flux_.access();
  }

  const XT::Common::ParameterType& parameter_type() const override final
  {
    return numerical_flux_.access().parameter_type();
  }

  void apply(const DF& source, DF& range, const XT::Common::Parameter& mu = {}) const
  {
    range.vector() *= 0.;
    LocalizableOperatorBase<GL, DF, DF> walker(grid_layer_, source, range);
    walker.append(local_coupling_operator_, new XT::Grid::ApplyOn::InnerIntersectionsPrimally<GL>(), mu);
    walker.append(local_coupling_operator_, new XT::Grid::ApplyOn::PeriodicIntersectionsPrimally<GL>(), mu);
    walker.walk();
  }

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
  const LocalCouplingOperatorType local_coupling_operator_;
}; // class AdvectionFvOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_FV_HH
