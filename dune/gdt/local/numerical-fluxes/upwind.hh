// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH

#include "interface.hh"

namespace Dune {
namespace GDT {


template <size_t d, size_t m = 1, class R = double>
class NumericalUpwindFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>
{
public:
  template <class... Args>
  explicit NumericalUpwindFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>()
  {}
};

template <size_t d, class R>
class NumericalUpwindFlux<d, 1, R> : public NumericalFluxInterface<d, 1, R>
{
  static const constexpr size_t m = 1;
  using ThisType = NumericalUpwindFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  NumericalUpwindFlux(const FluxType& flx)
    : BaseType(flx)
  {}

  NumericalUpwindFlux(const ThisType& other) = default;

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
    const auto df = this->flux().jacobian((u + v) / 2., param);
    if ((n * df) > 0)
      return this->flux().evaluate(u, param) * n;
    else
      return this->flux().evaluate(v, param) * n;
  }
}; // class NumericalUpwindFlux


template <size_t d, size_t m, class R>
NumericalUpwindFlux<d, m, R> make_numerical_upwind_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalUpwindFlux<d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH
