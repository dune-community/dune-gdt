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


template <size_t d, size_t m = 1, class R = double>
class NumericalEngquistOsherFlux : public internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>
{
public:
  template <class... Args>
  explicit NumericalEngquistOsherFlux(Args&&... /*args*/)
    : internal::ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>()
  {}
};

template <size_t d, class R>
class NumericalEngquistOsherFlux<d, 1, R> : public NumericalFluxInterface<d, 1, R>
{
  static const constexpr size_t m = 1;
  using ThisType = NumericalEngquistOsherFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  NumericalEngquistOsherFlux(const FluxType& flx)
    : BaseType(flx)
  {}

  NumericalEngquistOsherFlux(const ThisType& other) = default;

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
    auto integrate_f = [&](const auto& s, const std::function<double(const R&, const R&)>& min_max) {
      if (!(s[0] > 0.))
        return 0.;
      double ret = 0.;
      const OneDGrid state_grid(1, 0., s[0]);
      const auto state_interval = *state_grid.leafGridView().template begin<0>();
      for (const auto& quadrature_point :
           QuadratureRules<R, 1>::rule(state_interval.type(), this->flux().order(param))) {
        const auto local_uu = quadrature_point.position();
        const auto uu = state_interval.geometry().global(local_uu);
        const auto df = this->flux().jacobian(uu, param);
        ret += state_interval.geometry().integrationElement(local_uu) * quadrature_point.weight() * min_max(n * df, 0.);
      }
      return ret;
    };
    return (this->flux().evaluate(0., param) * n)
           + integrate_f(u, [](const double& a, const double& b) { return std::max(a, b); })
           + integrate_f(v, [](const double& a, const double& b) { return std::min(a, b); });
  }
}; // class NumericalEngquistOsherFlux


template <size_t d, size_t m, class R>
NumericalEngquistOsherFlux<d, m, R>
make_numerical_engquist_osher_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux)
{
  return NumericalEngquistOsherFlux<d, m, R>(flux);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_ENGQUIST_OSHER_HH
