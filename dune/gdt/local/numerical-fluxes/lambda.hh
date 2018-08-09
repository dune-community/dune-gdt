// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_numerical_flux_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_numerical_flux_HH

#include <functional>

#include "interface.hh"

namespace Dune {
namespace GDT {


/**
 * \brief Implementation of NumericalFluxInterface for a given lambda expression.
 */
template <size_t d, size_t m = 1, class R = double>
class NumericalLambdaFlux : public NumericalFluxInterface<d, m, R>
{
  using ThisType = NumericalLambdaFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  using LambdaType = std::function<StateRangeType(
      const StateRangeType&, const StateRangeType&, const PhysicalDomainType&, const XT::Common::Parameter&)>;

  NumericalLambdaFlux(const FluxType& flx, LambdaType lambda, const XT::Common::ParameterType& param_type = {})
    : BaseType(flx, param_type)
    , numerical_flux_(lambda)
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  StateRangeType apply(const StateRangeType& u,
                       const StateRangeType& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& param = {}) const override final
  {
    return numerical_flux_(u, v, n, this->parse_parameter(param));
  }

private:
  const LambdaType numerical_flux_;
}; // class NumericalLambdaFlux


template <size_t d, size_t m, class R>
NumericalLambdaFlux<d, m, R> make_numerical_lambda_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux,
                                                        typename NumericalLambdaFlux<d, m, R>::LambdaType lambda,
                                                        const XT::Common::ParameterType& param_type = {})
{
  return NumericalLambdaFlux<d, m, R>(flux, lambda, param_type);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_numerical_flux_HH
