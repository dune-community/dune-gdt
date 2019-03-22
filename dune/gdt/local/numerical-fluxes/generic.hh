// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_GENERIC_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_GENERIC_HH

#include <functional>

#include "interface.hh"

namespace Dune {
namespace GDT {


/**
 * \brief Implementation of NumericalFluxInterface for a given lambda expression.
 */
template <class I, size_t d, size_t m = 1, class R = double>
class GenericNumericalFlux : public NumericalFluxInterface<I, d, m, R>
{
  using ThisType = GenericNumericalFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::FunctionType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;

  using GenericFunctionType = std::function<StateType(const Intersection&,
                                                      const LocalIntersectionCoords&,
                                                      const StateType&,
                                                      const StateType&,
                                                      const PhysicalDomainType&,
                                                      const XT::Common::Parameter&)>;

  GenericNumericalFlux(const FluxType& flx, GenericFunctionType func, const XT::Common::ParameterType& param_type = {})
    : BaseType(flx, param_type)
    , numerical_flux_(func)
  {}

  GenericNumericalFlux(const FunctionType& flx,
                       GenericFunctionType func,
                       const XT::Common::ParameterType& param_type = {})
    : BaseType(flx, param_type)
    , numerical_flux_(func)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  StateType apply(const I& intersection,
                  const LocalIntersectionCoords& x,
                  const StateType& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const override final
  {
    return numerical_flux_(intersection, x, u, v, n, this->parse_parameter(param));
  }

private:
  const GenericFunctionType numerical_flux_;
}; // class GenericNumericalFlux


template <class I, size_t d, size_t m, class R>
GenericNumericalFlux<I, d, m, R>
make_generic_numerical_flux(const XT::Functions::FluxFunctionInterface<I, m, d, m, R>& flux,
                            typename GenericNumericalFlux<I, d, m, R>::GenericFunctionType func,
                            const XT::Common::ParameterType& param_type = {})
{
  return GenericNumericalFlux<I, d, m, R>(flux, func, param_type);
}

template <class I, size_t d, size_t m, class R>
GenericNumericalFlux<I, d, m, R>
make_generic_numerical_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux,
                            typename GenericNumericalFlux<I, d, m, R>::GenericFunctionType func,
                            const XT::Common::ParameterType& param_type = {})
{
  return GenericNumericalFlux<I, d, m, R>(flux, func, param_type);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_GENERIC_HH
