// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_ENGQUIST_OSHER_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_ENGQUIST_OSHER_HH

#include <functional>

#include <dune/geometry/quadraturerules.hh>

#include <dune/grid/onedgrid.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class I, size_t d, size_t m = 1, class R = double>
class NumericalEngquistOsherFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>
{
public:
  template <class... Args>
  explicit NumericalEngquistOsherFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<I, d, m, R>()
  {}
};

template <class I, size_t d, class R>
class NumericalEngquistOsherFlux<I, d, 1, R> : public NumericalFluxInterface<I, d, 1, R>
{
  static const constexpr size_t m = 1;
  using ThisType = NumericalEngquistOsherFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::LocalFluxType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;
  using typename BaseType::XIndependentFluxType;

  NumericalEngquistOsherFlux(const FluxType& flx)
    : BaseType(flx)
  {}

  NumericalEngquistOsherFlux(const XIndependentFluxType& flx)
    : BaseType(flx)
  {}

  NumericalEngquistOsherFlux(const ThisType& other) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  StateType apply(const LocalIntersectionCoords& x_in_intersection_coords,
                  const StateType& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const override final
  {
    this->compute_entity_coords(x_in_intersection_coords);
    auto integrate_f = [&](const LocalFluxType& local_flux,
                           const auto& x,
                           const auto& s,
                           const std::function<double(const R&, const R&)>& min_max) {
      if (!(s[0] > 0.))
        return 0.;
      double ret = 0.;
      const OneDGrid state_grid(1, 0., s[0]);
      const auto state_interval = *state_grid.leafGridView().template begin<0>();
      for (const auto& quadrature_point : QuadratureRules<R, 1>::rule(state_interval.type(), local_flux.order(param))) {
        const auto local_uu = quadrature_point.position();
        const auto uu = state_interval.geometry().global(local_uu);
        const auto df = local_flux.jacobian(x, uu, param);
        ret += state_interval.geometry().integrationElement(local_uu) * quadrature_point.weight() * min_max(n * df, 0.);
      }
      return ret;
    };
    return (local_flux_inside_->evaluate(x_in_inside_coords_, 0., param) * n)
           + integrate_f(*local_flux_inside_,
                         x_in_inside_coords_,
                         u,
                         [](const double& a, const double& b) { return std::max(a, b); })
           + integrate_f(this->intersection().neighbor() ? *local_flux_outside_ : *local_flux_inside_,
                         x_in_outside_coords_,
                         v,
                         [](const double& a, const double& b) { return std::min(a, b); });
  }

private:
  using BaseType::local_flux_inside_;
  using BaseType::local_flux_outside_;
  using BaseType::x_in_inside_coords_;
  using BaseType::x_in_outside_coords_;
}; // class NumericalEngquistOsherFlux


template <class I, size_t d, size_t m, class R>
NumericalEngquistOsherFlux<I, d, m, R>
make_numerical_engquist_osher_flux(const XT::Functions::FluxFunctionInterface<I, m, d, m, R>& flux)
{
  return NumericalEngquistOsherFlux<I, d, m, R>(flux);
}

template <class I, size_t d, size_t m, class R>
NumericalEngquistOsherFlux<I, d, m, R>
make_numerical_engquist_osher_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalEngquistOsherFlux<I, d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_ENGQUIST_OSHER_HH
