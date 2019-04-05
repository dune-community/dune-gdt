// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH

#include <cmath>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class I, size_t d, size_t m = 1, class R = double>
class NumericalLaxFriedrichsFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>
{
public:
  template <class... Args>
  explicit NumericalLaxFriedrichsFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>()
  {}
};

template <class I, size_t d, class R>
class NumericalLaxFriedrichsFlux<I, d, 1, R> : public NumericalFluxInterface<I, d, 1, R>
{
  static const constexpr size_t m = 1;
  using ThisType = NumericalLaxFriedrichsFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;
  using typename BaseType::XIndependentFluxType;

  NumericalLaxFriedrichsFlux(const FluxType& flx)
    : BaseType(flx)
  {}

  NumericalLaxFriedrichsFlux(const XIndependentFluxType& func)
    : BaseType(func)
  {}

  NumericalLaxFriedrichsFlux(const ThisType& other) = default;

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
    // prepare
    mutable_this->bind(intersection);
    this->compute_entity_coords(intersection, x);
    // evaluate
    const auto df_u = local_flux_inside_->jacobian(x_in_inside_coords_, u, param);
    const auto df_v = local_flux_outside_->jacobian(x_in_outside_coords_, v, param);
    const auto f_u = local_flux_inside_->evaluate(x_in_inside_coords_, u, param);
    const auto f_v = local_flux_outside_->evaluate(x_in_outside_coords_, v, param);
    // compute numerical flux
    const auto lambda = 1. / std::max(df_u.infinity_norm(), df_v.infinity_norm());
    return 0.5 * ((f_u + f_v) * n + (u - v) / lambda);
  }

private:
  using BaseType::local_flux_inside_;
  using BaseType::local_flux_outside_;
  using BaseType::mutable_this;
  using BaseType::x_in_inside_coords_;
  using BaseType::x_in_outside_coords_;
}; // class NumericalLaxFriedrichsFlux


template <class I, size_t d, size_t m, class R>
NumericalLaxFriedrichsFlux<I, d, m, R>
make_numerical_lax_friedrichs_flux(const XT::Functions::FluxFunctionInterface<I, m, d, m, R>& flux)
{
  return NumericalLaxFriedrichsFlux<I, d, m, R>(flux);
}

template <class I, size_t d, size_t m, class R>
NumericalLaxFriedrichsFlux<I, d, m, R>
make_numerical_lax_friedrichs_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalLaxFriedrichsFlux<I, d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_LAX_FRIEDRICHS_HH
