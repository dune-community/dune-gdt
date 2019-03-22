// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/flux-function.hh>
#include <dune/xt/functions/interfaces/function.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


/**
 * Given the sought solution of a system of m conservation laws, u: R^d -> R^m, and d flux
 * functions f_s: R^d x R^m -> R^m, for 1 <= s <= d (modelled by the
 * flux f: R^d x R^m -> R^{d x m}), the purpose of a numerical flux
 * g: R^m x R^m x R^d -> R^m is to approximate f(.) * n, e.g., g(x, u, u, n) = f(x, u) * n.
 */
template <class Intersection, size_t d, size_t m = 1, class R = double>
class NumericalFluxInterface : public XT::Common::ParametricInterface
{
  using ThisType = NumericalFluxInterface;

public:
  using I = Intersection;
  using E = typename I::Entity;
  using FluxType = XT::Functions::FluxFunctionInterface<E, m, d, m, R>;
  using FunctionType = XT::Functions::FunctionInterface<m, d, m, R>;
  using FunctionWrapperType = XT::Functions::StateFunctionAsFluxFunctionWrapper<E, m, d, m, R>;
  using PhysicalDomainType = typename FluxType::DomainType;
  using LocalIntersectionCoords = FieldVector<typename I::ctype, d - 1>;
  using StateType = typename FluxType::StateType;

  NumericalFluxInterface(const FluxType& flx, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx.parameter_type())
    , flux_(flx)
    , x_dependent_(true)
  {}

  NumericalFluxInterface(FluxType*&& flx_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx_ptr->parameter_type())
    , flux_(flx_ptr)
    , x_dependent_(true)
  {}

  NumericalFluxInterface(const FunctionType& func, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + func.parameter_type())
    , flux_(new FunctionWrapperType(func))
    , x_dependent_(false)
  {}

  NumericalFluxInterface(FunctionType*&& func_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + func_ptr->parameter_type())
    , flux_(new FunctionWrapperType(func_ptr))
    , x_dependent_(false)
  {}


  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual bool linear() const
  {
    return false;
  }

  virtual bool x_dependent() const
  {
    return x_dependent_;
  }

  const FluxType& flux() const
  {
    return flux_.access();
  }

  void compute_entity_coords(const I& intersection, const LocalIntersectionCoords& x_in_local_intersection_coords) const
  {
    if (this->x_dependent()) {
      x_in_inside_coords_ =
          intersection.inside().geometry().local(intersection.geometry().global(x_in_local_intersection_coords));
      if (intersection.neighbor())
        x_in_outside_coords_ =
            intersection.outside().geometry().local(intersection.geometry().global(x_in_local_intersection_coords));
      else
        x_in_outside_coords_ = x_in_inside_coords_;
    }
  }

  virtual StateType apply(const I& intersection,
                          const LocalIntersectionCoords& x_in_local_intersection_coords,
                          const StateType& u,
                          const StateType& v,
                          const PhysicalDomainType& n,
                          const XT::Common::Parameter& param = {}) const = 0;

  template <class V>
  StateType apply(const I& intersection,
                  const LocalIntersectionCoords x_in_local_intersection_coords,
                  const StateType& u,
                  const XT::LA::VectorInterface<V>& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(v.size() != m, Exceptions::numerical_flux_error, "v.size() = " << v.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      v_[ii] = v[ii];
    return this->apply(intersection, x_in_local_intersection_coords, u, v_, n, param);
  }

  template <class U>
  StateType apply(const I& intersection,
                  const LocalIntersectionCoords x_in_local_intersection_coords,
                  const XT::LA::VectorInterface<U>& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(u.size() != m, Exceptions::numerical_flux_error, "u.size() = " << u.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      u_[ii] = u[ii];
    return this->apply(intersection, x_in_local_intersection_coords, u_, v, n, param);
  }

  template <class U, class V>
  StateType apply(const I& intersection,
                  const LocalIntersectionCoords x_in_local_intersection_coords,
                  const XT::LA::VectorInterface<U>& u,
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
    return this->apply(intersection, x_in_local_intersection_coords, u_, v_, n, param);
  } // ... apply(...)

private:
  const XT::Common::ConstStorageProvider<FluxType> flux_;
  const bool x_dependent_;
  mutable StateType u_;
  mutable StateType v_;

protected:
  mutable PhysicalDomainType x_in_inside_coords_;
  mutable PhysicalDomainType x_in_outside_coords_;
}; // class NumericalFluxInterface


namespace internal {


template <class Intersection, size_t d, size_t m = 1, class R = double>
class ThisNumericalFluxIsNotAvailableForTheseDimensions : public NumericalFluxInterface<Intersection, d, m, R>
{
  using ThisType = ThisNumericalFluxIsNotAvailableForTheseDimensions;
  using BaseType = NumericalFluxInterface<Intersection, d, m, R>;

public:
  using typename BaseType::I;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;

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

  StateType apply(const I& /*intersection*/,
                  const PhysicalDomainType /*x_in_local_intersection_coords*/,
                  const StateType& /*u*/,
                  const StateType& /*v*/,
                  const PhysicalDomainType& /*n*/,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW(Exceptions::numerical_flux_error, "d = " << d << "\n   m = " << m);
    return StateType();
  }
}; // class ThisNumericalFluxIsNotAvailableForTheseDimensions


} // namespace internal


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_INTERFACE_HH
