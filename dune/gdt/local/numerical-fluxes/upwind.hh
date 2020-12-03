// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class I, size_t d, size_t m = 1, class R = double>
class NumericalUpwindFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>
{
public:
  template <class... Args>
  explicit NumericalUpwindFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>()
  {}
};

template <class I, size_t d, class R>
class NumericalUpwindFlux<I, d, 1, R> : public NumericalFluxInterface<I, d, 1, R>
{
  static constexpr size_t m = 1;
  using ThisType = NumericalUpwindFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;
  using typename BaseType::XIndependentFluxType;

  NumericalUpwindFlux(const FluxType& flx)
    : BaseType(flx)
  {}

  NumericalUpwindFlux(const XIndependentFluxType& func)
    : BaseType(func)
  {}

  NumericalUpwindFlux(const ThisType& other) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  StateType apply(const LocalIntersectionCoords& x,
                  const StateType& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const override final
  {
    this->compute_entity_coords(x);
    const auto df = local_flux_inside_->jacobian(x_in_inside_coords_, (u + v) / 2., param);
    if (n * df[0] > 0)
      return local_flux_inside_->evaluate(x_in_inside_coords_, u, param) * n;
    else
      return local_flux_outside_->evaluate(x_in_outside_coords_, v, param) * n;
  }

private:
  using BaseType::local_flux_inside_;
  using BaseType::local_flux_outside_;
  using BaseType::x_in_inside_coords_;
  using BaseType::x_in_outside_coords_;
}; // class NumericalUpwindFlux


template <class E, size_t d, size_t m, class R>
auto make_numerical_upwind_flux(const XT::Functions::FluxFunctionInterface<E, m, d, m, R>& flux)
{
  using I = XT::Grid::extract_entity_t<E>;
  return NumericalUpwindFlux<I, d, m, R>(flux);
}

template <class I, // <- has to be specified manually
          size_t d,
          size_t m,
          class R>
auto make_numerical_upwind_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalUpwindFlux<I, d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_UPWIND_HH
