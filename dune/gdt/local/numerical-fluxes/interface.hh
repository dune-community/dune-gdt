// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/function.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


/**
 * Given the sought solution of a system of m conservation laws, u: R^d -> R^m, and d flux functions f_s: R^m -> R^m,
 * for 1 <= s <= d (modelled by the flux f: R^m -> R^{d x m}), the purpose of a numerical flux
 * g: R^m x R^m x R^d -> R^m is to approximate f(.) * n, e.g., g(u, u, n) = f(u) * n.
 */
template <size_t d, size_t m = 1, class R = double>
class NumericalFluxInterface : public XT::Common::ParametricInterface
{
  using ThisType = NumericalFluxInterface<d, m, R>;

public:
  using FluxType = XT::Functions::FunctionInterface<m, d, m, R>;
  using PhysicalDomainType = FieldVector<double, d>;
  using StateRangeType = typename FluxType::DomainType;

  NumericalFluxInterface(const FluxType& flx, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx.parameter_type())
    , flux_(flx)
  {}

  NumericalFluxInterface(FluxType*&& flx_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx_ptr->parameter_type())
    , flux_(flx_ptr)
  {}

  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual bool linear() const
  {
    return false;
  }

  const FluxType& flux() const
  {
    return flux_.access();
  }

  virtual StateRangeType apply(const StateRangeType& u,
                               const StateRangeType& v,
                               const PhysicalDomainType& n,
                               const XT::Common::Parameter& param = {}) const = 0;

  template <class V>
  StateRangeType apply(const StateRangeType& u,
                       const XT::LA::VectorInterface<V>& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(v.size() != m, Exceptions::numerical_flux_error, "v.size() = " << v.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      v_[ii] = v[ii];
    return this->apply(u, v_, n, param);
  }

  template <class U>
  StateRangeType apply(const XT::LA::VectorInterface<U>& u,
                       const StateRangeType& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(u.size() != m, Exceptions::numerical_flux_error, "u.size() = " << u.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      u_[ii] = u[ii];
    return this->apply(u_, v, n, param);
  }

  template <class U, class V>
  StateRangeType apply(const XT::LA::VectorInterface<U>& u,
                       const XT::LA::VectorInterface<V>& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(u.size() != m, Exceptions::numerical_flux_error, "u.size() = " << u.size() << "\n   m = " << m);
    DUNE_THROW_IF(v.size() != m, Exceptions::numerical_flux_error, "v.size() = " << v.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii) {
      u_[ii] = u[ii];
      v_[ii] = v[ii];
    }
    return this->apply(u_, v_, n, param);
  } // ... apply(...)

private:
  const XT::Common::ConstStorageProvider<FluxType> flux_;
  mutable StateRangeType u_;
  mutable StateRangeType v_;
}; // class NumericalFluxInterface


namespace internal {


template <size_t d, size_t m = 1, class R = double>
class ThisNumericalFluxIsNotAvailableForTheseDimensions : public NumericalFluxInterface<d, m, R>
{
  using ThisType = ThisNumericalFluxIsNotAvailableForTheseDimensions<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  template <class... Args>
  explicit ThisNumericalFluxIsNotAvailableForTheseDimensions(Args&&... /*args*/)
    : BaseType(new XT::Functions::ConstantFunction<m, d, m, R>(0.))
  {
    DUNE_THROW(Exceptions::numerical_flux_error, "d = " << d << "\n   m = " << m);
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    DUNE_THROW(Exceptions::numerical_flux_error, "d = " << d << "\n   m = " << m);
    return nullptr;
  }

  bool linear() const override final
  {
    DUNE_THROW(Exceptions::numerical_flux_error, "d = " << d << "\n   m = " << m);
    return false;
  }

  using BaseType::apply;

  StateRangeType apply(const StateRangeType& /*u*/,
                       const StateRangeType& /*v*/,
                       const PhysicalDomainType& /*n*/,
                       const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW(Exceptions::numerical_flux_error, "d = " << d << "\n   m = " << m);
    return StateRangeType();
  }
}; // class ThisNumericalFluxIsNotAvailableForTheseDimensions


} // namespace internal


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH
