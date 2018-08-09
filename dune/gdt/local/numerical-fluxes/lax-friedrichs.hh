// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH

#include <cmath>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <size_t d, size_t m = 1, class R = double>
class NumericalLaxFriedrichsFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>
{
public:
  template <class... Args>
  explicit NumericalLaxFriedrichsFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>()
  {
  }
};

template <size_t d, class R>
class NumericalLaxFriedrichsFlux<d, 1, R> : public NumericalFluxInterface<d, 1, R>
{
  static const constexpr size_t m = 1;
  using ThisType = NumericalLaxFriedrichsFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  NumericalLaxFriedrichsFlux(const FluxType& flx)
    : BaseType(flx)
  {
  }

  NumericalLaxFriedrichsFlux(const ThisType& other) = default;

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
    const auto lambda =
        1. / std::max(this->flux().jacobian(u, param).infinity_norm(), this->flux().jacobian(v, param).infinity_norm());
    return 0.5 * ((this->flux().evaluate(u, param) + this->flux().evaluate(v, param)) * n) + 0.5 * ((u - v) / lambda);
  }
}; // class NumericalLaxFriedrichsFlux


template <size_t d, size_t m, class R>
NumericalLaxFriedrichsFlux<d, m, R>
make_numerical_lax_friedrichs_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalLaxFriedrichsFlux<d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH
