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
#include <dune/xt/common/vector.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>
#include <dune/xt/functions/interfaces/localizable-flux-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/tools/euler.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class SpaceType>
class LocalAdvectionFvCouplingOperator;

template <class SpaceType>
class LocalAdvectionFvBoundaryOperator;


namespace internal {


template <class SpaceType>
class LocalAdvectionFvCouplingOperatorTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgInnerOperator instead!");

public:
  using derived_type = LocalAdvectionFvCouplingOperator<SpaceType>;
};


template <class SpaceType>
class LocalAdvectionFvBoundaryOperatorTraits
{
  static_assert(is_fv_space<SpaceType>::value, "Use LocalAdvectionDgCouplingOperator instead!");

public:
  using derived_type = LocalAdvectionFvBoundaryOperator<SpaceType>;
};


} // namespace internal


template <class E, class D, size_t d, class R, size_t m>
class NumericalFluxInterface : public XT::Common::ParametricInterface
{
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
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalUpwindingFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

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
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalLaxFriedrichsFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

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
  using BaseType = NumericalFluxInterface<E, D, d, R, 1>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalEngquistOsherFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

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
 * \todo \attention Does not work for d > 1, the computation of df * n and the eigenvalue decomposition is probably
 *                  broken!
 * \note Checks can be disabled (to improve performance) by defining
 *       DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS
 */
template <class E, class D, size_t d, class R, size_t m>
class NumericalVijayasundaramFlux : public NumericalFluxInterface<E, D, d, R, m>
{
  static_assert(d == 1, "Unreliable for other dimension!");
  using BaseType = NumericalFluxInterface<E, D, d, R, m>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalVijayasundaramFlux(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    check_values(u);
    check_values(v);
    check_values(n);
    // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
    const auto df = this->flux().partial_u({}, 0.5 * (u + v));
    const auto P = convert_to_eigen_matrix(df * n);
    check_values(P);
    FieldVector<R, m> result(0.);
    try {
      // compute decomposition
      auto eigensolver = XT::LA::make_eigen_solver(P);
      const auto evs = eigensolver.real_eigenvalues();
      check_values(evs);
      const auto T = eigensolver.real_eigenvectors_as_matrix();
      check_values(T);
      auto T_inv = T.copy();
      try {
        T_inv = XT::LA::invert_matrix(T);
        check_values(T_inv);
      } catch (XT::LA::Exceptions::matrix_invert_failed_bc_result_is_not_a_left_inverse& ee) {
        const auto identity = XT::LA::eye_matrix<XT::LA::EigenDenseMatrix<R>>(m, m);
        // check each eigenvalue/eigenvector pair individually
        for (size_t ii = 0; ii < m; ++ii) {
          const auto& lambda = evs[ii];
          std::cerr << "\n\nchecking lambda_" << ii << " = " << lambda << ":" << std::endl;
          // try all possible eigenvectors
          XT::LA::EigenDenseVector<R> eigenvector(m, 0.);
          for (size_t jj = 0; jj < m; ++jj) {
            // get jth column of T
            for (size_t kk = 0; kk < m; ++kk)
              eigenvector[kk] = T.get_entry(kk, jj);
            std::cerr << "  testing column " << jj << " = " << eigenvector << ": ";
            const auto tolerance = 1e-15;
            const auto error = ((P - identity * lambda) * eigenvector).sup_norm();
            if (error < tolerance)
              std::cerr << "this IS an eigenvector (up to " << tolerance << ")!" << std::endl;
            else {
              std::cerr << "this IS NOT an eigenvector (|| (P - lambda I) * ev ||_infty = " << error << ")!"
                        << std::endl;
            }
          }
        }
        DUNE_THROW(
            InvalidStateException,
            "Eigen decomposition of flux jacobians P not successfull, because T is not invertible (see below)!\n\n"
                << "P = "
                << P
                << "\n\ndiagonal of lambda (eigenvalues) = "
                << evs
                << "\n\nT (eigenvectors) = "
                << T
                << "\n\n\n\nThis was the original error: "
                << ee.what());
      }

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS
      // test decomposition
      XT::LA::EigenDenseMatrix<R> lambda(m, m);
      for (size_t ii = 0; ii < m; ++ii)
        lambda.set_entry(ii, ii, evs[ii]);
      if (((T * lambda * T_inv) - P).sup_norm() / P.sup_norm() > 1e-14) {
        const auto identity = XT::LA::eye_matrix<XT::LA::EigenDenseMatrix<R>>(m, m);
        // check each eigenvalue/eigenvector pair individually
        for (size_t ii = 0; ii < m; ++ii) {
          const auto& lmbd = evs[ii];
          std::cerr << "\n\nchecking lambda_" << ii << " = " << lmbd << ":" << std::endl;
          // try all possible eigenvectors
          XT::LA::EigenDenseVector<R> eigenvector(m, 0.);
          for (size_t jj = 0; jj < m; ++jj) {
            // get jth column of T
            for (size_t kk = 0; kk < m; ++kk)
              eigenvector[kk] = T.get_entry(kk, jj);
            std::cerr << "  testing column " << jj << " = " << eigenvector << ": ";
            const auto tolerance = 1e-15;
            const auto error = ((P - identity * lmbd) * eigenvector).sup_norm();
            if (error < tolerance)
              std::cerr << "this IS an eigenvector (up to " << tolerance << ")!" << std::endl;
            else {
              std::cerr << "this IS NOT an eigenvector (|| (P - lambda I) * ev ||_infty = " << error << ")!"
                        << std::endl;
            }
          }
        }
        DUNE_THROW(InvalidStateException,
                   "Eigen decomposition of flux jacobians P not successfull!\n\n"
                       << "P = "
                       << P
                       << "\n\ndiagonal of lambda (eigenvalues) = "
                       << evs
                       << "\n\nT (eigenvectors) = "
                       << T
                       << "\n\n||(T * lambda * T_inv) - P||_infty / ||P||_infty = "
                       << ((T * lambda * T_inv) - P).sup_norm() / P.sup_norm());
      }
#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS

      // compute numerical flux [DF2016, p. 428, (8.108)]
      XT::LA::EigenDenseMatrix<R> lambda_plus(m, m, 0.);
      XT::LA::EigenDenseMatrix<R> lambda_minus(m, m, 0.);
      for (size_t ii = 0; ii < m; ++ii) {
        lambda_plus.set_entry(ii, ii, std::max(evs[ii], 0.));
        lambda_minus.set_entry(ii, ii, std::min(evs[ii], 0.));
      }
      check_values(lambda_plus);
      check_values(lambda_minus);
      const auto P_plus = convert_to_common_matrix(T * lambda_plus * T_inv);
      const auto P_minus = convert_to_common_matrix(T * lambda_minus * T_inv);
      check_values(P_plus);
      check_values(P_minus);

      result = P_plus * u + P_minus * v;
      check_values(result);
    } catch (XT::LA::Exceptions::eigen_solver_failed& ee) {
      DUNE_THROW(InvalidStateException,
                 "Eigen decomposition of flux jacobians P not successfull (see below for the original error)!\n\n"
                     << "u = "
                     << u
                     << "\n\nv = "
                     << v
                     << "\n\nDF((u + v)/2) = "
                     << df
                     << "\n\nP = "
                     << P
                     << "\n\nThis was the original error:\n\n"
                     << ee.what());
    }
    return result;
  } // ... apply(...)

private:
  template <class V>
  static typename std::enable_if<XT::Common::is_vector<V>::value, void>::type check_values(const V& vec)
  {
    for (size_t ii = 0; ii < vec.size(); ++ii)
      if (XT::Common::isnan(vec[ii]) || XT::Common::isinf(vec[ii]))
        DUNE_THROW(InvalidStateException, vec);
  }

  template <class M>
  static typename std::enable_if<XT::Common::is_matrix<M>::value, void>::type check_values(const M& mat)
  {
    using MM = XT::Common::MatrixAbstraction<M>;
    for (size_t ii = 0; ii < MM::rows(mat); ++ii)
      for (size_t jj = 0; jj < MM::cols(mat); ++jj)
        if (XT::Common::isnan(MM::get_entry(mat, ii, jj)) || XT::Common::isinf(MM::get_entry(mat, ii, jj)))
          DUNE_THROW(InvalidStateException, mat);
  }

  static XT::LA::EigenDenseMatrix<R> convert_to_eigen_matrix(const XT::Common::FieldMatrix<R, m, m>& source)
  {
    using Abstraction = XT::Common::MatrixAbstraction<XT::Common::FieldMatrix<R, m, m>>;
    XT::LA::EigenDenseMatrix<R> target(Abstraction::rows(source), Abstraction::cols(source), 0.);
    for (size_t ii = 0; ii < target.rows(); ++ii)
      for (size_t jj = 0; jj < target.cols(); ++jj)
        target.set_entry(ii, jj, Abstraction::get_entry(source, ii, jj));
    return target;
  }

  static XT::Common::FieldMatrix<R, m, m> convert_to_common_matrix(const XT::LA::EigenDenseMatrix<R>& source)
  {
    using TargetAbstraction = XT::Common::MatrixAbstraction<XT::Common::FieldMatrix<R, m, m>>;
    using SourceAbstraction = XT::Common::MatrixAbstraction<XT::LA::EigenDenseMatrix<R>>;
    auto target = TargetAbstraction::create(SourceAbstraction::rows(source), SourceAbstraction::cols(source), 0.);
    for (size_t ii = 0; ii < SourceAbstraction::rows(source); ++ii)
      for (size_t jj = 0; jj < SourceAbstraction::cols(source); ++jj)
        TargetAbstraction::set_entry(target, ii, jj, SourceAbstraction::get_entry(source, ii, jj));
    return target;
  }

  /// This was the first attempt to use an eigenvalue decomposition, but the one from eigen is bad!
  /// In addition, the P comin from df * n is also bad!
  //  auto vijayasundaram_euler = [&](const auto& u_, const auto& v_, const auto& n, const auto& /*mu*/) {
  //    if (d != 2)
  //      DUNE_THROW(NotImplemented, "Only for 2d!\nd = " << d);
  //    check_values(u_);
  //    check_values(v_);
  //    check_values(n);
  //    // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
  //    const auto w_conservative = 0.5 * (u_ + v_);
  //    const auto df = flux.partial_u({}, w_conservative);
  //    const auto P_generic = convert_to<XT::LA::EigenDenseMatrix<R>>(df * n);
  //    check_values(P_generic);

  //    // compute P matrix directly for euler [DF2016, p.405, (8.20)]
  //    const auto w_primitiv = to_primitive(w_conservative);
  //    const auto& rho = w_conservative[0];
  //    DomainType v;
  //    v[0] = w_primitiv[1];
  //    v[1] = w_primitiv[2];
  //    const auto& e = w_conservative[3];
  //    const auto gamma_1 = gamma - 1;
  //    const auto gamma_2 = gamma - 2;
  //    const auto G = gamma * e / rho - 0.5 * gamma_1 * v.two_norm2();
  //    XT::LA::EigenDenseMatrix<R> P_euler(m, m, 0.);
  //    P_euler.set_entry(0, 0, 0);
  //    P_euler.set_entry(0, 1, n[0]);
  //    P_euler.set_entry(0, 2, n[1]);
  //    P_euler.set_entry(0, 3, 0);

  //    P_euler.set_entry(1, 0, 0.5 * gamma_1 * v.two_norm2() * n[0] - v[0] * (v * n));
  //    P_euler.set_entry(1, 1, -gamma_2 * v[0] * n[0] + v * n);
  //    P_euler.set_entry(1, 2, v[0] * n[1] - gamma_1 * v[1] * n[0]);
  //    P_euler.set_entry(1, 3, gamma_1 * n[0]);

  //    P_euler.set_entry(2, 0, 0.5 * gamma_1 * v.two_norm2() * n[1] - v[1] * (v * n));
  //    P_euler.set_entry(2, 1, v[1] * n[0] - gamma_1 * v[0] * n[1]);
  //    P_euler.set_entry(2, 2, -gamma_2 * v[1] * n[1] + v * n);
  //    P_euler.set_entry(2, 3, gamma_1 * n[1]);

  //    P_euler.set_entry(3, 0, (gamma_1 * v.two_norm2() - gamma * e / rho) * (v * n));
  //    P_euler.set_entry(3, 1, G * n[0] - gamma_1 * v[0] * (v * n));
  //    P_euler.set_entry(3, 2, G * n[1] - gamma_1 * v[1] * (v * n));
  //    P_euler.set_entry(3, 3, gamma * (v * n));

  //    if (XT::Common::FloatCmp::ne(P_generic, P_euler))
  //      DUNE_THROW(InvalidStateException,
  //                 "Jacobian flux matrices do not coincide!\n\nThis is the generic P:\n\n"
  //                     << P_generic
  //                     << "\n\nThis is the 2d euler P:\n\n"
  //                     << P_euler
  //                     << "\n\nu = "
  //                     << u_
  //                     << "\n\nv = "
  //                     << v_
  //                     << "\n\nn = "
  //                     << n);

  //    const auto& P = P_euler;

  //    FieldVector<R, m> result(0.);
  //    try {
  //      // compute decomposition
  //      auto eigensolver = XT::LA::make_eigen_solver(P);
  //      const auto evs = eigensolver.real_eigenvalues();
  //      check_values(evs);
  //      const auto T = eigensolver.real_eigenvectors_as_matrix();
  //      check_values(T);
  //      auto T_inv = T.copy();
  //      try {
  //        T_inv = XT::LA::invert_matrix(T);
  //        //        check_values(T_inv);
  //      } catch (XT::LA::Exceptions::matrix_invert_failed& ee) {
  //        const auto identity = XT::LA::eye_matrix<XT::LA::EigenDenseMatrix<R>>(m, m);
  //        // check each eigenvalue/eigenvector pair individually
  //        for (size_t ii = 0; ii < m; ++ii) {
  //          const auto& lambda = evs[ii];
  //          std::cerr << "\n\nchecking lambda_" << ii << " = " << lambda << ":" << std::endl;
  //          // try all possible eigenvectors
  //          XT::LA::EigenDenseVector<R> eigenvector(m, 0.);
  //          for (size_t jj = 0; jj < m; ++jj) {
  //            // get jth column of T
  //            for (size_t kk = 0; kk < m; ++kk)
  //              eigenvector[kk] = T.get_entry(kk, jj);
  //            std::cerr << "  testing column " << jj << " = " << eigenvector << ": ";
  //            const auto tolerance = 1e-15;
  //            const auto error = ((P - identity * lambda) * eigenvector).sup_norm();
  //            if (error < tolerance)
  //              std::cerr << "this IS an eigenvector (up to " << tolerance << ")!" << std::endl;
  //            else {
  //              std::cerr << "this IS NOT an eigenvector (|| (P - lambda I) * ev ||_infty = " << error << ")!"
  //                        << std::endl;
  //            }
  //          }
  //        }
  //        DUNE_THROW(
  //            InvalidStateException,
  //            "Eigen decomposition of flux jacobians P not successfull, because T is not invertible (see below)!\n\n"
  //                << "P = "
  //                << P
  //                << "\n\ndiagonal of lambda (eigenvalues) = "
  //                << evs
  //                << "\n\nT (eigenvectors) = "
  //                << T
  //                << "\n\n\n\nThis was the original error: "
  //                << ee.what());
  //      }

  //#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS
  //      // test decomposition
  //      XT::LA::EigenDenseMatrix<R> lambda(m, m);
  //      for (size_t ii = 0; ii < m; ++ii)
  //        lambda.set_entry(ii, ii, evs[ii]);
  //      if (((T * lambda * T_inv) - P).sup_norm() / P.sup_norm() > 1e-10) {
  //        const auto identity = XT::LA::eye_matrix<XT::LA::EigenDenseMatrix<R>>(m, m);
  //        // check each eigenvalue/eigenvector pair individually
  //        for (size_t ii = 0; ii < m; ++ii) {
  //          const auto& lmbd = evs[ii];
  //          std::cerr << "\n\nchecking lambda_" << ii << " = " << lmbd << ":" << std::endl;
  //          // try all possible eigenvectors
  //          XT::LA::EigenDenseVector<R> eigenvector(m, 0.);
  //          for (size_t jj = 0; jj < m; ++jj) {
  //            // get jth column of T
  //            for (size_t kk = 0; kk < m; ++kk)
  //              eigenvector[kk] = T.get_entry(kk, jj);
  //            std::cerr << "  testing column " << jj << " = " << eigenvector << ": ";
  //            const auto tolerance = 1e-15;
  //            const auto error = ((P - identity * lmbd) * eigenvector).sup_norm();
  //            if (error < tolerance)
  //              std::cerr << "this IS an eigenvector (up to " << tolerance << ")!" << std::endl;
  //            else {
  //              std::cerr << "this IS NOT an eigenvector (|| (P - lambda I) * ev ||_infty = " << error << ")!"
  //                        << std::endl;
  //            }
  //          }
  //        }
  //        DUNE_THROW(InvalidStateException,
  //                   "Eigen decomposition of flux jacobians P not successfull!\n\n"
  //                       << "P = "
  //                       << P
  //                       << "\n\ndiagonal of lambda (eigenvalues) = "
  //                       << evs
  //                       << "\n\nT (eigenvectors) = "
  //                       << T
  //                       << "\n\n||(T * lambda * T_inv) - P||_infty / ||P||_infty = "
  //                       << ((T * lambda * T_inv) - P).sup_norm() / P.sup_norm());
  //      }
  //#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS

  //      // compute numerical flux [DF2016, p. 428, (8.108)]
  //      XT::LA::EigenDenseMatrix<R> lambda_plus(m, m, 0.);
  //      XT::LA::EigenDenseMatrix<R> lambda_minus(m, m, 0.);
  //      for (size_t ii = 0; ii < m; ++ii) {
  //        lambda_plus.set_entry(ii, ii, std::max(evs[ii], 0.));
  //        lambda_minus.set_entry(ii, ii, std::min(evs[ii], 0.));
  //      }
  //      check_values(lambda_plus);
  //      check_values(lambda_minus);
  //      const auto P_plus = convert_to<XT::Common::FieldMatrix<R, m, m>>(T * lambda_plus * T_inv);
  //      const auto P_minus = convert_to<XT::Common::FieldMatrix<R, m, m>>(T * lambda_minus * T_inv);
  //      check_values(P_plus);
  //      check_values(P_minus);

  //      result = P_plus * u_ + P_minus * v_;
  //      check_values(result);
  //    } catch (XT::LA::Exceptions::eigen_solver_failed& ee) {
  //      DUNE_THROW(InvalidStateException,
  //                 "Eigen decomposition of flux jacobians P not successfull (see below for the original error)!\n\n"
  //                     << "u = "
  //                     << u_
  //                     << "\n\nv = "
  //                     << v_
  //                     << "\n\nDF((u + v)/2) = "
  //                     << df
  //                     << "\n\nP = "
  //                     << P
  //                     << "\n\nThis was the original error:\n\n"
  //                     << ee.what());
  //    }
  //    return result;
  //  }; // ... vijayasundaram_euler(...)
}; // class NumericalVijayasundaramFlux


template <class E, class D, size_t d, class R, size_t m>
NumericalVijayasundaramFlux<E, D, d, R, m> make_numerical_vijayasundaram_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux)
{
  return NumericalVijayasundaramFlux<E, D, d, R, m>(flux);
}


template <class E, class D, size_t d, class R, size_t m>
class NumericalVijayasundaramEulerFlux
{
  static_assert(AlwaysFalse<E>::value, "Not implemented for these dimensions yet!");
};


/**
 * \note Checks can be disabled (to improve performance) by defining
 *       DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS
 */
template <class E, class D, class R>
class NumericalVijayasundaramEulerFlux<E, D, 2, R, 4> : public NumericalFluxInterface<E, D, 2, R, 4>
{
  static const constexpr size_t d = 2;
  static const constexpr size_t m = d + 2;
  using BaseType = NumericalFluxInterface<E, D, d, R, m>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit NumericalVijayasundaramEulerFlux(const double& gamma,
                                            const double& eigenvalue_check_tolerance,
                                            Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , euler_tools_(gamma)
    , tolerance_(eigenvalue_check_tolerance)
  {
  }

  RangeType apply(const RangeType& u,
                  const RangeType& v,
                  const DomainType& n,
                  const XT::Common::Parameter& /*mu*/ = {}) const override final
  {
    const auto w_conservative = 0.5 * (u + v);

    const auto eigenvalues = euler_tools_.eigenvalues_flux_jacobi_matrix(w_conservative, n);
    const auto eigenvectors = euler_tools_.eigenvectors_flux_jacobi_matrix(w_conservative, n);
    const auto eigenvectors_inv = euler_tools_.eigenvectors_inv_flux_jacobi_matrix(w_conservative, n);

#ifndef DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS
    const auto identity = XT::LA::eye_matrix<XT::Common::FieldMatrix<R, m, m>>(m, m);
    if ((eigenvectors_inv * eigenvectors - identity).infinity_norm() > tolerance_)
      DUNE_THROW(InvalidStateException,
                 "\n\neigenvectors:\n\n"
                     << eigenvectors
                     << "\n\neigenvectors_inverse:\n\n"
                     << eigenvectors_inv
                     << "\n\n|| eigenvectors_inv * eigenvectors - identity ||_infty = "
                     << (eigenvectors_inv * eigenvectors - identity).infinity_norm());

    const auto eigenvaluematrix = euler_tools_.eigenvaluematrix_flux_jacobi_matrix(w_conservative, n);
    if (((eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w_conservative, n) * eigenvectors)) - eigenvaluematrix)
            .infinity_norm()
        > tolerance_)
      DUNE_THROW(InvalidStateException,
                 "\n\neigenvectors:\n\n"
                     << eigenvectors
                     << "\n\neigenvectors_inverse:\n\n"
                     << eigenvectors_inv
                     << "\n\neigenvalues:"
                     << eigenvalues
                     << "\n\nP:\n\n"
                     << euler_tools_.flux_jacobi_matrix(w_conservative, n)
                     << "\n\neigenvectors_inv * (P * eigenvectors):\n\n"
                     << eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w_conservative, n) * eigenvectors)
                     << "\n\n|| eigenvectors_inv * (P * eigenvectors) - eigenvalues||_infty = "
                     << ((eigenvectors_inv * (euler_tools_.flux_jacobi_matrix(w_conservative, n) * eigenvectors))
                         - eigenvaluematrix)
                            .infinity_norm());
#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_DISABLE_CHECKS

    XT::Common::FieldMatrix<R, m, m> eigenvalues_plus(0.);
    for (size_t ii = 0; ii < m; ++ii)
      eigenvalues_plus[ii][ii] = std::max(0., eigenvalues[ii]);

    XT::Common::FieldMatrix<R, m, m> eigenvalues_minus(0.);
    for (size_t ii = 0; ii < m; ++ii)
      eigenvalues_minus[ii][ii] = std::min(0., eigenvalues[ii]);

    const auto P_plus = eigenvectors * eigenvalues_plus * eigenvectors_inv;
    const auto P_minus = eigenvectors * eigenvalues_minus * eigenvectors_inv;
    return P_plus * u + P_minus * v;
  } // ... apply(...)

private:
  const EulerTools<d, R> euler_tools_;
  const double tolerance_;
}; // class NumericalVijayasundaramEulerFlux


template <class E, class D, size_t d, class R, size_t m>
NumericalVijayasundaramEulerFlux<E, D, d, R, m> make_numerical_vijayasundaram_euler_flux(
    const XT::Functions::
        GlobalFluxFunctionInterface<E, D, d, XT::Functions::LocalizableFunctionInterface<E, D, d, R, m, 1>, 0, R, d, m>&
            flux,
    const double& gamma,
    const double eigenvalue_check_tolerance = 1e-10)
{
  return NumericalVijayasundaramEulerFlux<E, D, d, R, m>(gamma, eigenvalue_check_tolerance, flux);
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
      DUNE_THROW(InvalidStateException, "");
    if (v.size() != u_outside->vector().size())
      DUNE_THROW(InvalidStateException, "");
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


template <class SpaceType>
class LocalAdvectionFvBoundaryTreatmentInterface : public XT::Common::ParametricInterface
{
  static_assert(is_space<SpaceType>::value, "");
  static_assert(SpaceType::dimRangeCols == 1, "Not Implemented yet!");
  using E = typename SpaceType::EntityType;
  using D = typename SpaceType::DomainFieldType;
  static const constexpr size_t d = SpaceType::dimDomain;
  using R = typename SpaceType::RangeFieldType;
  static const constexpr size_t r = SpaceType::dimRange;
  static const constexpr size_t rC = 1;
  using StateType = XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>;

public:
  using IntersectionType = XT::Grid::extract_intersection_t<typename SpaceType::GridLayerType>;
  using IntersectionPointType = FieldVector<D, d - 1>;
  using RangeType = typename StateType::RangeType;
  using NumericalFluxType = NumericalFluxInterface<E, D, d, R, r>;
  using FluxType = typename NumericalFluxType::FluxType;

  LocalAdvectionFvBoundaryTreatmentInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalAdvectionFvBoundaryTreatmentInterface() = default;

  virtual RangeType apply(const IntersectionType& intersection,
                          const IntersectionPointType& x_intersection,
                          const FluxType& flux,
                          const RangeType& u,
                          const XT::Common::Parameter& mu = {}) const = 0;
}; // class LocalAdvectionFvBoundaryTreatmentInterface


template <class SpaceType>
class LocalAdvectionFvLambdaBoundaryTreatment : public LocalAdvectionFvBoundaryTreatmentInterface<SpaceType>
{
  using BaseType = LocalAdvectionFvBoundaryTreatmentInterface<SpaceType>;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::IntersectionPointType;
  using typename BaseType::RangeType;
  using typename BaseType::NumericalFluxType;
  using typename BaseType::FluxType;

  using LambdaType = std::function<RangeType(const IntersectionType&,
                                             const IntersectionPointType&,
                                             const FluxType&,
                                             const RangeType&,
                                             const XT::Common::Parameter&)>;

  LocalAdvectionFvLambdaBoundaryTreatment(LambdaType lambda, const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , lambda_(lambda)
  {
  }

  RangeType apply(const IntersectionType& intersection,
                  const IntersectionPointType& x_intersection,
                  const FluxType& flux,
                  const RangeType& u,
                  const XT::Common::Parameter& mu = {}) const override final
  {
    return lambda_(intersection, x_intersection, flux, u, this->parse_parameter(mu));
  }

private:
  const LambdaType lambda_;
}; // class LocalAdvectionFvLambdaBoundaryTreatment


/**
 * \note Presumes that the basis evaluates to 1.
 * \todo Improve local vector conversion in apply.
 */
template <class SpaceType>
class LocalAdvectionFvBoundaryOperator
    : public LocalBoundaryOperatorInterface<internal::LocalAdvectionFvBoundaryOperatorTraits<SpaceType>>
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
  using BoundaryTreatmentType = LocalAdvectionFvBoundaryTreatmentInterface<SpaceType>;

  LocalAdvectionFvBoundaryOperator(const NumericalFluxType& numerical_flux,
                                   const BoundaryTreatmentType& boundary_treatment)
    : numerical_flux_(numerical_flux)
    , boundary_treatment_(boundary_treatment)
  {
  }

  template <class VectorType, class I>
  void apply(const ConstDiscreteFunction<SpaceType, VectorType>& source,
             const I& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range) const
  {
    static_assert(XT::Grid::is_intersection<I>::value, "");
    const auto& entity = local_range.entity();
    const auto u_inside = source.local_discrete_function(entity);
    const auto x_intersection = ReferenceElements<D, d - 1>::general(intersection.type()).position(0, 0);
    const auto normal = intersection.unitOuterNormal(x_intersection);
    // copy the local DoF vector to matching FieldVector
    typename StateType::RangeType u;
    if (u.size() != u_inside->vector().size())
      DUNE_THROW(InvalidStateException, "");
    for (size_t ii = 0; ii < u.size(); ++ii)
      u[ii] = u_inside->vector().get(ii);
    const auto v = boundary_treatment_.apply(intersection, x_intersection, numerical_flux_.flux(), u, /*mu*/ {});
    const auto g = numerical_flux_.apply(u, v, normal, /*mu*/ {});
    const auto h = entity.geometry().volume();
    for (size_t ii = 0; ii < r; ++ii)
      local_range.vector().add(ii, (g[ii] * intersection.geometry().volume()) / h);
  } // ... apply(...)

private:
  const NumericalFluxType& numerical_flux_;
  const BoundaryTreatmentType& boundary_treatment_;
}; // class LocalAdvectionFvBoundaryOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_ADVECTION_FV_HH
