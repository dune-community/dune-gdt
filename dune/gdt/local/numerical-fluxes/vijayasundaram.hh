// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_VIJAYASUNDARAM_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_VIJAYASUNDARAM_HH

#include <functional>
#include <tuple>

#include <dune/xt/common/matrix.hh>
#include <dune/xt/common/math.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


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
                               XT::Common::FieldMatrix<XT::Common::real_t<R>, m, m>>(
          const FluxType&,
          const FieldVector<R, m>&,
          const FieldVector<double, d>&,
          const XT::Common::Parameter& param)>;

  NumericalVijayasundaramFlux(const FluxType& flx)
    : BaseType(flx)
    , flux_eigen_decomposition_([](const auto& f, const auto& w, const auto& n, const auto& param) {
      // evaluate flux jacobian, compute P matrix [DF2016, p. 404, (8.17)]
      const auto df = f.jacobian(w, param);
      const auto P = df * n;
      auto eigensolver = XT::LA::make_eigen_solver(
          P, {{"type", XT::LA::eigen_solver_types(P)[0]}, {"assert_real_eigendecomposition", "1e-10"}});
      return std::make_tuple(
          eigensolver.real_eigenvalues(), eigensolver.real_eigenvectors(), eigensolver.real_eigenvectors_inverse());
    })
  {
  }

  NumericalVijayasundaramFlux(const FluxType& flx, FluxEigenDecompositionLambdaType flux_eigen_decomposition)
    : BaseType(flx)
    , flux_eigen_decomposition_(flux_eigen_decomposition)
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
    // compute decomposition
    const auto eigendecomposition = flux_eigen_decomposition_(this->flux(), 0.5 * (u + v), n, param);
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
  const FluxEigenDecompositionLambdaType flux_eigen_decomposition_;
}; // class NumericalVijayasundaramFlux


template <size_t d, size_t m, class R, class... Args>
NumericalVijayasundaramFlux<d, m, R>
make_numerical_vijayasundaram_flux(const XT::Functions::FunctionInterface<m, d, m, R>& flux, Args&&... args)
{
  return NumericalVijayasundaramFlux<d, m, R>(flux, std::forward<Args>(args)...);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_VIJAYASUNDARAM_HH
