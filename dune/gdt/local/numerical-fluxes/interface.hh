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
class NumericalFluxInterface
  : public XT::Grid::IntersectionBoundObject<Intersection>
  , public XT::Common::ParametricInterface
{
  using ThisType = NumericalFluxInterface;

public:
  using I = Intersection;
  using E = typename I::Entity;
  using FluxType = XT::Functions::FluxFunctionInterface<E, m, d, m, R>;
  using LocalFluxType = typename FluxType::LocalFunctionType;
  using XIndependentFluxType = XT::Functions::FunctionInterface<m, d, m, R>;
  using FluxWrapperType = XT::Functions::StateFunctionAsFluxFunctionWrapper<E, m, d, m, R>;
  using PhysicalDomainType = typename FluxType::DomainType;
  using LocalIntersectionCoords = FieldVector<typename I::ctype, d - 1>;
  using StateType = typename FluxType::StateType;
  using DynamicStateType = typename FluxType::DynamicStateType;

  NumericalFluxInterface(const FluxType& flx, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx.parameter_type())
    , flux_(flx)
    , local_flux_inside_(flux_.access().local_function())
    , local_flux_outside_(flux_.access().local_function())
  {}

  NumericalFluxInterface(FluxType*&& flx_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx_ptr->parameter_type())
    , flux_(flx_ptr)
    , local_flux_inside_(flux_.access().local_function())
    , local_flux_outside_(flux_.access().local_function())
  {}

  NumericalFluxInterface(const XIndependentFluxType& func, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + func.parameter_type())
    , flux_(new FluxWrapperType(func))
    , local_flux_inside_(flux_.access().local_function())
    , local_flux_outside_(flux_.access().local_function())
  {}

  NumericalFluxInterface(XIndependentFluxType*&& func_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + func_ptr->parameter_type())
    , flux_(new FluxWrapperType(func_ptr))
    , local_flux_inside_(flux_.access().local_function())
    , local_flux_outside_(flux_.access().local_function())
  {}

  NumericalFluxInterface(const ThisType& other)
    : XT::Grid::IntersectionBoundObject<Intersection>(other)
    , XT::Common::ParametricInterface(other)
    , flux_(other.flux_)
    , local_flux_inside_(flux_.access().local_function())
    , local_flux_outside_(flux_.access().local_function())
  {}

  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual bool linear() const
  {
    return false;
  }

  virtual bool x_dependent() const
  {
    return flux_.access().x_dependent();
  }

  const FluxType& flux() const
  {
    return flux_.access();
  }

  // One of the two following apply methods has to be implemented by derived classes
  virtual StateType apply(const LocalIntersectionCoords& x_in_local_intersection_coords,
                          const StateType& u,
                          const StateType& v,
                          const PhysicalDomainType& n,
                          const XT::Common::Parameter& param = {}) const
  {
    DynamicStateType ret(m, 0.);
    apply(x_in_local_intersection_coords, u, v, n, ret, param);
    return XT::Common::convert_to<StateType>(ret);
  }

  virtual void apply(const LocalIntersectionCoords& x_in_local_intersection_coords,
                     const DynamicStateType& u,
                     const DynamicStateType& v,
                     const PhysicalDomainType& n,
                     DynamicStateType& ret,
                     const XT::Common::Parameter& param = {}) const
  {
    ret = XT::Common::convert_to<DynamicStateType>(apply(x_in_local_intersection_coords, u, v, n, param));
  }

  // Convenience apply methods
  template <class V>
  StateType apply(const LocalIntersectionCoords x_in_local_intersection_coords,
                  const StateType& u,
                  const XT::LA::VectorInterface<V>& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(v.size() != m, Exceptions::numerical_flux_error, "v.size() = " << v.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      v_[ii] = v[ii];
    return this->apply(x_in_local_intersection_coords, u, v_, n, param);
  }

  template <class U>
  StateType apply(const LocalIntersectionCoords x_in_local_intersection_coords,
                  const XT::LA::VectorInterface<U>& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(u.size() != m, Exceptions::numerical_flux_error, "u.size() = " << u.size() << "\n   m = " << m);
    for (size_t ii = 0; ii < m; ++ii)
      u_[ii] = u[ii];
    return this->apply(x_in_local_intersection_coords, u_, v, n, param);
  }

  template <class U, class V>
  StateType apply(const LocalIntersectionCoords x_in_local_intersection_coords,
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
    return this->apply(x_in_local_intersection_coords, u_, v_, n, param);
  } // ... apply(...)

  using XT::Grid::IntersectionBoundObject<I>::intersection;

private:
  const XT::Common::ConstStorageProvider<FluxType> flux_;
  mutable StateType u_;
  mutable StateType v_;

protected:
  void post_bind(const I& inter) override
  {
    local_flux_inside_->bind(inter.inside());
    if (inter.neighbor())
      local_flux_outside_->bind(inter.outside());
  }

  void compute_entity_coords(const LocalIntersectionCoords& x_in_local_intersection_coords) const
  {
    if (this->x_dependent()) {
      if (!this->is_bound_)
        DUNE_THROW(Dune::InvalidStateException, "You have to call bind(intersection) before calling this function!");
      x_in_inside_coords_ =
          intersection().inside().geometry().local(intersection().geometry().global(x_in_local_intersection_coords));
      if (intersection().neighbor())
        x_in_outside_coords_ =
            intersection().outside().geometry().local(intersection().geometry().global(x_in_local_intersection_coords));
      else
        x_in_outside_coords_ = x_in_inside_coords_;
    }
  }

  mutable std::unique_ptr<LocalFluxType> local_flux_inside_;
  mutable std::unique_ptr<LocalFluxType> local_flux_outside_;
  mutable PhysicalDomainType x_in_inside_coords_;
  mutable PhysicalDomainType x_in_outside_coords_;
}; // class NumericalFluxInterface


namespace internal {


template <class Intersection, size_t d, size_t m = 1, class R = double>
class ThisNumericalFluxIsNotAvailableForTheseDimensions : public NumericalFluxInterface<Intersection, d, m, R>
{
  using BaseType = NumericalFluxInterface<Intersection, d, m, R>;

public:
  using typename BaseType::I;
  using typename BaseType::LocalIntersectionCoords;
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

  StateType apply(const LocalIntersectionCoords& /*x_in_local_intersection_coords*/,
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
