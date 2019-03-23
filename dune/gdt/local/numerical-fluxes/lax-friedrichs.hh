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


template <class I, size_t d, size_t m, class R = double>
class NumericalLaxFriedrichsFlux : public NumericalFluxInterface<I, d, m, R>
{
  using ThisType = NumericalLaxFriedrichsFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;
  using typename BaseType::XIndependentFluxType;
  using FluxJacobianType = XT::Common::FieldVector<XT::Common::FieldMatrix<R, m, m>, d>;

  NumericalLaxFriedrichsFlux(const FluxType& flx, const double lambda = 0.)
    : BaseType(flx)
    , lambda_(lambda)
  {
    if (XT::Common::is_zero(lambda_) && m != 1)
      DUNE_THROW(Dune::NotImplemented, "Not yet implemented for m > 1 if lambda is not provided!");
  }

  NumericalLaxFriedrichsFlux(const XIndependentFluxType& func, const double lambda = 0.)
    : BaseType(func)
    , lambda_(lambda)
  {
    if (XT::Common::is_zero(lambda_) && m != 1)
      DUNE_THROW(Dune::NotImplemented, "Not yet implemented for m > 1 if lambda is not provided!");
  }

  NumericalLaxFriedrichsFlux(const ThisType& other) = default;

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
    // prepare
    this->compute_entity_coords(x);
    // evaluate
    R lambda = lambda_;
    if (XT::Common::is_zero(lambda)) {
      const auto df_u = local_flux_inside_->jacobian(x_in_inside_coords_, u, param);
      const auto df_v = local_flux_outside_->jacobian(x_in_outside_coords_, v, param);
      for (size_t dd = 0; dd < d; ++dd) {
        lambda = std::max(lambda, df_u[dd].infinity_norm());
        lambda = std::max(lambda, df_v[dd].infinity_norm());
      }
      lambda = 1./lambda;
    }
    const auto f_u = local_flux_inside_->evaluate(x_in_inside_coords_, u, param);
    const auto f_v = local_flux_outside_->evaluate(x_in_outside_coords_, v, param);
    StateType ret(0.);
    for (size_t dd = 0; dd < d; ++dd)
      ret += (f_u[dd] + f_v[dd]) * (n[dd] * 0.5);
    ret += (u - v) * (0.5 / lambda);
    return ret;
  }

private:
  using BaseType::local_flux_inside_;
  using BaseType::local_flux_outside_;
  using BaseType::x_in_inside_coords_;
  using BaseType::x_in_outside_coords_;
  using BaseType::mutable_this;
  const double lambda_;
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
