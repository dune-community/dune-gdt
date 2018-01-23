// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH

#include <functional>
#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/onedgrid.hh>

#include <dune/xt/common/densevector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/vector.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/matrix-inverter.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/tools/euler.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class SpaceType>
class LocalAdvectionFvCouplingOperator;

template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomExtrapolation;

template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux;


namespace internal {


template <class SpaceType>
class LocalAdvectionFvCouplingOperatorTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgInnerOperator instead!");

public:
  using derived_type = LocalAdvectionFvCouplingOperator<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomExtrapolationTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgCouplingOperatorByCustomExtrapolation instead!");

public:
  using derived_type = LocalAdvectionFvBoundaryOperatorByCustomExtrapolation<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgCouplingOperatorByCustomFlux instead!");

public:
  using derived_type = LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux<SpaceType>;
};


} // namespace internal


template <class E, class D, size_t d, class R, size_t m>
class NumericalFluxInterface : public XT::Common::ParametricInterface
{
  using ThisType = NumericalFluxInterface<E, D, d, R, m>;

public:
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>;
  using FluxType = XT::Functions::GlobalFluxFunctionInterface<E, D, d, StateType, 0, R, d, m>;
  using DomainType = typename StateType::DomainType;
  using RangeType = typename StateType::RangeType;

  NumericalFluxInterface(const FluxType& flx, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , flux_(flx)
  {
  }

  NumericalFluxInterface(FluxType*&& flx_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
    , flux_(flx_ptr)
  {
  }

  NumericalFluxInterface(const ThisType&) = default;
  NumericalFluxInterface(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  const FluxType& flux() const
  {
    return flux_.access();
  }

  virtual RangeType
  apply(const RangeType& u, const RangeType& v, const DomainType& n, const XT::Common::Parameter& mu = {}) const = 0;

private:
  const XT::Common::ConstStorageProvider<FluxType> flux_;
}; // class NumericalFluxInterface


template <class LF>
class NumericalLambdaFlux : public NumericalFluxInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, LF::r>
{
  static_assert(XT::Functions::is_localizable_function<LF>::value, "");
  static_assert(LF::rC == 1, "");
  using ThisType = NumericalLambdaFlux<LF>;
  using BaseType = NumericalFluxInterface<typename LF::E, typename LF::D, LF::d, typename LF::R, LF::r>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  using LambdaType =
      std::function<RangeType(const RangeType&, const RangeType&, const DomainType&, const XT::Common::Parameter&)>;

  NumericalLambdaFlux(const FluxType& flx, LambdaType lambda, const XT::Common::ParameterType& param_type = {})
    : BaseType(flx, param_type)
    , lambda_(lambda)
  {
  }

  NumericalLambdaFlux(const ThisType&) = default;
  NumericalLambdaFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& mu = {}) const override final
  {
    return lambda_(u, v, n, this->parse_parameter(mu));
  }

private:
  const LambdaType lambda_;
}; // class NumericalLambdaFlux


template <class E, class D, size_t d, class R, size_t m>
class NumericalUpwindingFlux
{
  static_assert(AlwaysFalse<E>::value, "Not implemented for systems yet!");
};


template <class E, class D, size_t d, class R>
class NumericalUpwindingFlux<E, D, d, R, 1> : public NumericalFluxInterface<E, D, d, R, 1>
{
  using ThisType = NumericalUpwindingFlux<E, D, d, R, 1>;
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalUpwindingFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  NumericalUpwindingFlux(const ThisType&) = default;
  NumericalUpwindingFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    const auto df = this->flux().partial_u({}, (u + v) / 2.);
    if ((n * df) > 0)
      return this->flux().evaluate({}, u) * n;
    else
      return this->flux().evaluate({}, v) * n;
  }
}; // class NumericalUpwindingFlux


template <class E, class D, size_t d, class R, size_t m>
NumericalUpwindingFlux<E, D, d, R, m> make_numerical_upwinding_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux)
{
  return NumericalUpwindingFlux<E, D, d, R, m>(flux);
}


template <class E, class D, size_t d, class R, size_t m>
class NumericalLaxFriedrichsFlux
{
  static_assert(AlwaysFalse<E>::value, "Not implemented for systems yet!");
};


template <class E, class D, size_t d, class R>
class NumericalLaxFriedrichsFlux<E, D, d, R, 1> : public NumericalFluxInterface<E, D, d, R, 1>
{
  using ThisType = NumericalLaxFriedrichsFlux<E, D, d, R, 1>;
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalLaxFriedrichsFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  NumericalLaxFriedrichsFlux(const ThisType&) = default;
  NumericalLaxFriedrichsFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    const auto lambda =
        1. / std::max(this->flux().partial_u({}, u).infinity_norm(), this->flux().partial_u({}, v).infinity_norm());
    return 0.5 * ((this->flux().evaluate({}, u) + this->flux().evaluate({}, v)) * n) + 0.5 * ((u - v) / lambda);
  }
}; // class NumericalLaxFriedrichsFlux


template <class E, class D, size_t d, class R, size_t m>
NumericalLaxFriedrichsFlux<E, D, d, R, m> make_numerical_lax_friedrichs_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux)
{
  return NumericalLaxFriedrichsFlux<E, D, d, R, m>(flux);
}


template <class E, class D, size_t d, class R, size_t m>
class NumericalEngquistOsherFlux
{
  static_assert(AlwaysFalse<E>::value, "Not implemented for systems yet!");
};


template <class E, class D, size_t d, class R>
class NumericalEngquistOsherFlux<E, D, d, R, 1> : public NumericalFluxInterface<E, D, d, R, 1>
{
  using ThisType = NumericalEngquistOsherFlux<E, D, d, R, 1>;
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalEngquistOsherFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  NumericalEngquistOsherFlux(const ThisType&) = default;
  NumericalEngquistOsherFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    auto integrate_f = [&](const auto& s, const std::function<double(const R&, const R&)>& min_max) {
      if (!(s[0] > 0.))
        return 0.;
      D ret = 0.;
      const OneDGrid state_grid(1, 0., s[0]);
      const auto state_interval = *state_grid.leafGridView().template begin<0>();
      for (const auto& quadrature_point : QuadratureRules<R, 1>::rule(state_interval.type(), this->flux().order())) {
        const auto local_uu = quadrature_point.position();
        const auto uu = state_interval.geometry().global(local_uu);
        const auto df = this->flux().partial_u({}, uu);
        ret += state_interval.geometry().integrationElement(local_uu) * quadrature_point.weight() * min_max(n * df, 0.);
      }
      return ret;
    };
    return (this->flux().evaluate({}, 0.) * n)
           + integrate_f(u, [](const double& a, const double& b) { return std::max(a, b); })
           + integrate_f(v, [](const double& a, const double& b) { return std::min(a, b); });
  }
}; // class NumericalEngquistOsherFlux


template <class E, class D, size_t d, class R, size_t m>
NumericalEngquistOsherFlux<E, D, d, R, m> make_numerical_engquist_osher_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux)
{
  return NumericalEngquistOsherFlux<E, D, d, R, m>(flux);
}


/**
 * \note Checks can be disabled (to improve performance) by defining
 *       DUNE_GDT_DISABLE_CHECKS
 */
template <class E, class D, size_t d, class R, size_t m>
class NumericalVijayasundaramFlux : public NumericalFluxInterface<E, D, d, R, m>
{
  using ThisType = NumericalVijayasundaramFlux<E, D, d, R, m>;
  using BaseType = NumericalFluxInterface<E, D, d, R, m>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;
  using typename BaseType::FluxType;

  using FluxEigenDecompositionLambdaType =
      std::function<std::tuple<std::vector<XT::Common::real_t<R>>,
                               XT::Common::FieldMatrix<XT::Common::real_t<R>, m, m>,
                               XT::Common::FieldMatrix<XT::Common::real_t<R>, m, m>>(const FieldVector<R, m>&,
                                                                                     const FieldVector<D, d>&)>;

  NumericalVijayasundaramFlux(const FluxType& flx)
    : BaseType(flx)
    , flux_eigen_decomposition_lambda_([&](const auto& w, const auto& n) {
      // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
      const auto df = XT::Common::make_field_container(this->flux().partial_u({}, w));
      const auto P = df * n;
      auto eigensolver = XT::LA::make_eigen_solver(
          P, {{"type", XT::LA::eigen_solver_types(P)[0]}, {"assert_real_eigendecomposition", "1e-10"}});
      return std::make_tuple(
          eigensolver.real_eigenvalues(), eigensolver.real_eigenvectors(), eigensolver.real_eigenvectors_inverse());
    })
  {
  }

  NumericalVijayasundaramFlux(const FluxType& flx, FluxEigenDecompositionLambdaType flux_eigen_decomposition_lambda)
    : BaseType(flx)
    , flux_eigen_decomposition_lambda_(flux_eigen_decomposition_lambda)
  {
  }

  NumericalVijayasundaramFlux(const ThisType&) = default;
  NumericalVijayasundaramFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    // compute decomposition
    const auto eigendecomposition = flux_eigen_decomposition_lambda_(0.5 * (u + v), n);
    const auto& evs = std::get<0>(eigendecomposition);
    const auto& T = std::get<1>(eigendecomposition);
    const auto& T_inv = std::get<2>(eigendecomposition);
    // compute numerical flux [DF2016, p. 428, (8.108)]
    auto lambda_plus = XT::Common::zeros_like(T);
    auto lambda_minus = XT::Common::zeros_like(T);
    for (size_t ii = 0; ii < m; ++ii) {
      const auto& real_ev = evs[ii];
      XT::Common::set_matrix_entry(lambda_plus, ii, ii, std::max(real_ev, 0.));
      XT::Common::set_matrix_entry(lambda_minus, ii, ii, std::min(real_ev, 0.));
    }
    const auto P_plus = T * lambda_plus * T_inv;
    const auto P_minus = T * lambda_minus * T_inv;
    return P_plus * u + P_minus * v;
  } // ... apply(...)

private:
  const FluxEigenDecompositionLambdaType flux_eigen_decomposition_lambda_;
}; // class NumericalVijayasundaramFlux


template <class E, class D, size_t d, class R, size_t m, class... Args>
NumericalVijayasundaramFlux<E, D, d, R, m> make_numerical_vijayasundaram_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux,
    Args&&... args)
{
  return NumericalVijayasundaramFlux<E, D, d, R, m>(flux, std::forward<Args>(args)...);
}


/**
 * \note Checks can be disabled (to improve performance) by defining DUNE_GDT_DISABLE_CHECKS
 */
template <class E, class D, size_t d, class R>
class NumericalVijayasundaramEulerFlux : public NumericalVijayasundaramFlux<E, D, d, R, d + 2>
{
  static const constexpr size_t m = d + 2;
  using ThisType = NumericalVijayasundaramEulerFlux<E, D, d, R>;
  using BaseType = NumericalVijayasundaramFlux<E, D, d, R, m>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  explicit NumericalVijayasundaramEulerFlux(const FluxType& flx,
                                            const double& gamma,
                                            const double& eigenvalue_check_tolerance)
    : BaseType(flx,
               [&](const auto& w, const auto& n) {
                 const auto eigenvalues = euler_tools_.eigenvalues_flux_jacobi_matrix(w, n);
                 const auto eigenvectors = euler_tools_.eigenvectors_flux_jacobi_matrix(w, n);
                 const auto eigenvectors_inv = euler_tools_.eigenvectors_inv_flux_jacobi_matrix(w, n);

#ifndef DUNE_GDT_DISABLE_CHECKS
                 const auto identity = XT::LA::eye_matrix<XT::Common::FieldMatrix<R, m, m>>(m, m);
                 if ((eigenvectors_inv * eigenvectors - identity).infinity_norm() > tolerance_)
                   DUNE_THROW(Exceptions::numerical_flux_error,
                              "\n\neigenvectors:\n\n"
                                  << eigenvectors
                                  << "\n\neigenvectors_inverse:\n\n"
                                  << eigenvectors_inv
                                  << "\n\n|| eigenvectors_inv * eigenvectors - identity ||_infty = "
                                  << (eigenvectors_inv * eigenvectors - identity).infinity_norm());

                 const auto eigenvaluematrix = euler_tools_.eigenvaluematrix_flux_jacobi_matrix(w, n);
                 if (((eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w, n) * eigenvectors)) - eigenvaluematrix)
                         .infinity_norm()
                     > tolerance_)
                   DUNE_THROW(Exceptions::numerical_flux_error,
                              "\n\neigenvectors:\n\n"
                                  << eigenvectors
                                  << "\n\neigenvectors_inverse:\n\n"
                                  << eigenvectors_inv
                                  << "\n\neigenvalues:"
                                  << eigenvalues
                                  << "\n\nP:\n\n"
                                  << euler_tools_.flux_jacobi_matrix(w, n)
                                  << "\n\neigenvectors_inv * (P * eigenvectors):\n\n"
                                  << eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w, n) * eigenvectors)
                                  << "\n\n|| eigenvectors_inv * (P * eigenvectors) - eigenvalues||_infty = "
                                  << ((eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w, n) * eigenvectors))
                                      - eigenvaluematrix)
                                         .infinity_norm());
#endif // DUNE_GDT_DISABLE_CHECKS
                 return std::make_tuple(eigenvalues, eigenvectors, eigenvectors_inv);
               })
    , euler_tools_(gamma)
    , tolerance_(eigenvalue_check_tolerance)
  {
  }

  NumericalVijayasundaramEulerFlux(const ThisType&) = default;
  NumericalVijayasundaramEulerFlux(ThisType&&) = default;

  ThisType& operator=(const ThisType&) = delete;
  ThisType& operator=(ThisType&&) = delete;

private:
  const EulerTools<d, R> euler_tools_;
  const double tolerance_;
}; // class NumericalVijayasundaramEulerFlux


template <class E, class D, size_t d, class R>
NumericalVijayasundaramEulerFlux<E, D, d, R> make_numerical_vijayasundaram_euler_flux(
    const XT::Functions::GlobalFluxFunctionInterface<E,
                                                     D,
                                                     d,
                                                     XT::Functions::LocalizableFunctionInterface<E, D, d, R, d + 2, 1>,
                                                     0,
                                                     R,
                                                     d,
                                                     d + 2>& flux,
    const double& gamma,
    const double eigenvalue_check_tolerance = 1e-10)
{
  return NumericalVijayasundaramEulerFlux<E, D, d, R>(flux, gamma, eigenvalue_check_tolerance);
}


/**
 * \note Presumes that the basis evaluates to 1.
 * \todo Improve local vector handling in apply.
 */
template <class SpaceType>
class LocalAdvectionFvCouplingOperator
    : public LocalCouplingOperatorInterface<internal::LocalAdvectionFvCouplingOperatorTraits<SpaceType>>
{
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;

  LocalAdvectionFvCouplingOperator(const NumericalFluxType& numerical_flux)
    : numerical_flux_(numerical_flux)
  {
  }

  const XT::Common::ParameterType& parameter_type() const override final
  {
    return numerical_flux_.parameter_type();
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor,
             const XT::Common::Parameter& mu = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range_entity.entity();
    const auto& neighbor = local_range_neighbor.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto u_outside = source.local_discrete_function(neighbor);
    const auto normal = intersection.centerUnitOuterNormal();
    // copy the local DoF vector to matching FieldVectors
    typename StateType::RangeType u;
    typename StateType::RangeType v;
    if (u.size() != u_inside->vector().size())
      DUNE_THROW(Exceptions::local_operator_error,
                 "u.size() = " << u.size() << "\n   u_inside->vector().size() = " << u_inside->vector().size());
    if (v.size() != u_outside->vector().size())
      DUNE_THROW(Exceptions::local_operator_error,
                 "v.size() = " << v.size() << "\n   u_outside->vector().size() = " << u_outside->vector().size());
    for (size_t ii = 0; ii < u.size(); ++ii) {
      u[ii] = u_inside->vector().get(ii);
      v[ii] = u_outside->vector().get(ii);
    }
    const auto g = numerical_flux_.apply(u, v, normal, mu);
    const auto h = entity.geometry().volume();
    for (size_t ii = 0; ii < r; ++ii) {
      local_range_entity.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
      local_range_neighbor.vector().add(ii, (-1.0 * g[ii] * intersection.geometry().volume()) / h);
    }
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;
}; // class LocalAdvectionFvCouplingOperator


/**
 * \note Presumes that the basis evaluates to 1.
 * \todo Improve local vector conversion in apply.
 */
template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomExtrapolation
    : public LocalBoundaryOperatorInterface<internal::
                                                LocalAdvectionFvBoundaryOperatorByCustomExtrapolationTraits<SpaceType>>
{
  using BaseType =
      LocalBoundaryOperatorInterface<internal::LocalAdvectionFvBoundaryOperatorByCustomExtrapolationTraits<SpaceType>>;
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using IntersectionType = XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>;
  using IntersectionPointType = FieldVector<D, d - 1>;
  using RangeType = typename StateType::RangeType;
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;
  using FluxType = typename NumericalFluxType::FluxType;
  using LambdaType = std::function<RangeType(const IntersectionType&,
                                             const IntersectionPointType&,
                                             const FluxType&,
                                             const RangeType&,
                                             const XT::Common::Parameter&)>;

  LocalAdvectionFvBoundaryOperatorByCustomExtrapolation(
      const NumericalFluxType& numerical_flux,
      LambdaType boundary_treatment_lambda,
      const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(boundary_treatment_param_type)
    , numerical_flux_(numerical_flux)
    , boundary_treatment_(boundary_treatment_lambda)
  {
    if (numerical_flux_.parameter_type() != boundary_treatment_param_type)
      DUNE_THROW(NotImplemented, "todo: merge boundary_treatment_param_type and numerical_flux.parameter_type()!");
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range,
             const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto x_intersection = ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0);
    const auto normal = intersection.unitOuterNormal(x_intersection);
    // copy the local DoF vector to matching FieldVector
    typename StateType::RangeType u;
    if (u.size() != u_inside->vector().size())
      DUNE_THROW(Exceptions::local_operator_error,
                 "u.size() = " << u.size() << "\n   u_inside->vector().size() = " << u_inside->vector().size());
    for (size_t ii = 0; ii < u.size(); ++ii)
      u[ii] = u_inside->vector().get(ii);
    const auto v = boundary_treatment_(intersection, x_intersection, numerical_flux_.flux(), u, param);
    const auto g = numerical_flux_.apply(u, v, normal, param);
    const auto h = entity.geometry().volume();
    for (size_t ii = 0; ii < r; ++ii)
      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;
  const LambdaType boundary_treatment_;
}; // class LocalAdvectionFvBoundaryOperatorByCustomExtrapolation


/**
 * \note Presumes that the basis evaluates to 1.
 * \todo Improve local vector conversion in apply.
 */
template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux
    : public LocalBoundaryOperatorInterface<internal::
                                                LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxTraits<SpaceType>>
{
  using BaseType =
      LocalBoundaryOperatorInterface<internal::LocalAdvectionFvBoundaryOperatorByCustomNumericalFluxTraits<SpaceType>>;
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r>;

public:
  using DomainType = typename StateType::DomainType;
  using RangeType = typename StateType::RangeType;
  using LambdaType = std::function<RangeType(const RangeType&, const DomainType&, const XT::Common::Parameter&)>;

  LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux(
      LambdaType boundary_numerical_flux_lambda, const XT::Common::ParameterType& boundary_treatment_param_type = {})
    : BaseType(boundary_treatment_param_type)
    , boundary_numerical_flux_lambda_(boundary_numerical_flux_lambda)
  {
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range,
             const XT::Common::Parameter& param = {}) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto x_intersection = ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0);
    const auto normal = intersection.unitOuterNormal(x_intersection);
    // copy the local DoF vector to matching FieldVector
    typename StateType::RangeType u;
    if (u.size() != u_inside->vector().size())
      DUNE_THROW(Exceptions::local_operator_error,
                 "u.size() = " << u.size() << "\n   u_inside->vector().size() = " << u_inside->vector().size());
    for (size_t ii = 0; ii < u.size(); ++ii)
      u[ii] = u_inside->vector().get(ii);
    const auto g = boundary_numerical_flux_lambda_(u, normal, param);
    const auto h = entity.geometry().volume();
    for (size_t ii = 0; ii < r; ++ii)
      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  } // ... apply(...)

private:
  const LambdaType boundary_numerical_flux_lambda_;
}; // class LocalAdvectionFvBoundaryOperatorByCustomNumericalFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
