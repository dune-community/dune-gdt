// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH

#include <functional>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/functions/interfaces/function.hh>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"

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
  {
  }

  NumericalFluxInterface(FluxType*&& flx_ptr, const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type + flx_ptr->parameter_type())
    , flux_(flx_ptr)
  {
  }

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

  template <class U, class V>
  StateRangeType apply(const XT::LA::VectorInterface<U>& u,
                       const XT::LA::VectorInterface<V>& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(u.size() != m, Exceptions::operator_error, "u.size() = " << u.size() << "\n   m = " << m);
    DUNE_THROW_IF(v.size() != m, Exceptions::operator_error, "v.size() = " << v.size() << "\n   m = " << m);
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


/**
 * \brief Implementation of NumericalFluxInterface for a given lambda expression.
 */
template <size_t d, size_t m = 1, class R = double>
class NumericalLambdaFlux : public NumericalFluxInterface<d, m, R>
{
  using ThisType = NumericalLambdaFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  using LambdaType = std::function<StateRangeType(
      const StateRangeType&, const StateRangeType&, const PhysicalDomainType&, const XT::Common::Parameter&)>;

  NumericalLambdaFlux(const FluxType& flx, LambdaType lambda, const XT::Common::ParameterType& param_type = {})
    : BaseType(flx, param_type)
    , lambda_(lambda)
  {
  }

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
    return lambda_(u, v, n, this->parse_parameter(param));
  }

private:
  const LambdaType lambda_;
}; // class NumericalLambdaFlux


template <size_t d, size_t m, class R>
NumericalLambdaFlux<d, m, R> make_numerical_lambda_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux,
                                                        typename NumericalLambdaFlux<d, m, R>::LambdaType lambda,
                                                        const XT::Common::ParameterType& param_type = {})
{
  return NumericalLambdaFlux<d, m, R>(flux, lambda, param_type);
}


template <size_t d, size_t m = 1, class R = double>
class NumericalVijayasundaramFlux : public NumericalFluxInterface<d, m, R>
{
  using ThisType = NumericalVijayasundaramFlux<d, m, R>;
  using BaseType = NumericalFluxInterface<d, m, R>;

public:
  using typename BaseType::FluxType;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateRangeType;

  using FluxEigenDecompositionLambdaType =
      std::function<std::tuple<std::vector<XT::Common::real_t<R>>,
                               XT::Common::FieldMatrix<XT::Common::real_t<R>, m, m>,
                               XT::Common::FieldMatrix<XT::Common::real_t<R>, m, m>>(const FieldVector<R, m>&,
                                                                                     const FieldVector<double, d>&)>;

  NumericalVijayasundaramFlux(const FluxType& flx)
    : BaseType(flx)
    , flux_eigen_decomposition_lambda_([&](const auto& w, const auto& n) {
      // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
      const auto df = this->flux().jacobian(w);
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

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;

  StateRangeType apply(const StateRangeType& u,
                       const StateRangeType& v,
                       const PhysicalDomainType& n,
                       const XT::Common::Parameter& /*param*/ = {}) const override final
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
      XT::Common::set_matrix_entry(lambda_plus, ii, ii, XT::Common::max(real_ev, 0.));
      XT::Common::set_matrix_entry(lambda_minus, ii, ii, XT::Common::min(real_ev, 0.));
    }
    const auto P_plus = T * lambda_plus * T_inv;
    const auto P_minus = T * lambda_minus * T_inv;
    return P_plus * u + P_minus * v;
  } // ... apply(...)

private:
  FluxEigenDecompositionLambdaType flux_eigen_decomposition_lambda_;
}; // class NumericalVijayasundaramFlux


template <size_t d, size_t m, class R, class... Args>
NumericalVijayasundaramFlux<d, m, R>
make_numerical_vijayasundaram_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux, Args&&... args)
{
  return NumericalVijayasundaramFlux<d, m, R>(flux, std::forward<Args>(args)...);
}


/**
 * \note Presumes that the basis evaluates to 1.
 *
 * \note See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 */
template <class I,
          class SV,
          class SGV,
          size_t m = 1,
          class SR = double,
          class RR = SR,
          class IRGV = SGV,
          class IRV = SV,
          class ORR = RR,
          class ORGV = IRGV,
          class ORV = IRV>
class LocalAdvectionFvCouplingOperator
    : public LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>
{
  using ThisType = LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SR, RR, IRGV, IRV, ORR, ORGV, ORV>;
  using BaseType = LocalIntersectionOperatorInterface<I, SV, SGV, m, 1, SR, m, 1, RR, IRGV, IRV, ORGV, ORV>;

public:
  using BaseType::d;
  using typename BaseType::IntersectionType;
  using typename BaseType::SourceType;
  using typename BaseType::LocalInsideRangeType;
  using typename BaseType::LocalOutsideRangeType;

  using NumericalFluxType = NumericalFluxInterface<d, m, RR>;

  LocalAdvectionFvCouplingOperator(const NumericalFluxType& numerical_flux)
    : BaseType(numerical_flux.parameter_type())
    , numerical_flux_(numerical_flux.copy())
  {
  }

  LocalAdvectionFvCouplingOperator(const ThisType& other)
    : BaseType(other.parameter_type())
    , numerical_flux_(other.numerical_flux_->copy())
  {
  }

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalInsideRangeType& local_range_inside,
             LocalOutsideRangeType& local_range_outside,
             const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF((source.space().type() != SpaceType::finite_volume)
                      || (local_range_inside.space().type() != SpaceType::finite_volume)
                      || (local_range_outside.space().type() != SpaceType::finite_volume),
                  Exceptions::operator_error,
                  "Use LocalAdvectionDgCouplingOperator instead!");
    const auto& inside_element = local_range_inside.element();
    const auto& outside_element = local_range_outside.element();
    const auto u = source.local_discrete_function(inside_element);
    const auto v = source.local_discrete_function(outside_element);
    const auto normal = intersection.centerUnitOuterNormal();
    const auto g = numerical_flux_->apply(u->dofs(), v->dofs(), normal, param);
    const auto h_intersection = intersection.geometry().volume();
    for (size_t ii = 0; ii < m; ++ii) {
      local_range_inside.dofs()[ii] += (g[ii] * h_intersection) / inside_element.geometry().volume();
      local_range_outside.dofs()[ii] -= (g[ii] * h_intersection) / outside_element.geometry().volume();
    }
  } // ... apply(...)

private:
  std::unique_ptr<NumericalFluxType> numerical_flux_;
}; // class LocalAdvectionFvCouplingOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
