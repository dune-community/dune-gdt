// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2014, 2016 - 2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_SYMMETRIZED_LAPLACE_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_SYMMETRIZED_LAPLACE_HH

#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing scalar function lambda computes
 * `lambda(x) * 1/2*(\nabla phi(x) + (\nabla phi(x))^T) : \nabla psi(x)` for all combinations of phi in the ansatz and
 * psi in the test basis. Here, ':' denotes the (matrix) scalar product.
 */
template <class E, class F = double>
class LocalSymmetrizedLaplaceIntegrand
  : public LocalBinaryElementIntegrandInterface<E, E::dimension, 1, F, F, E::dimension, 1, F>
{

  using BaseType = LocalBinaryElementIntegrandInterface<E, E::dimension, 1, F, F, E::dimension, 1, F>;
  using ThisType = LocalSymmetrizedLaplaceIntegrand;

public:
  using BaseType::d;
  static constexpr size_t r = d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;

  LocalSymmetrizedLaplaceIntegrand(XT::Functions::GridFunction<E, 1, 1, F> diffusion_factor = F(1))
    : BaseType(diffusion_factor.parameter_type())
    , diffusion_factor_(diffusion_factor.copy_as_grid_function())
    , local_diffusion_factor_(diffusion_factor_->local_function())
  {}

  LocalSymmetrizedLaplaceIntegrand(const ThisType& other)
    : BaseType(other)
    , diffusion_factor_(other.diffusion_factor_->copy_as_grid_function())
    , local_diffusion_factor_(diffusion_factor_->local_function())
  {}

  LocalSymmetrizedLaplaceIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override
  {
    local_diffusion_factor_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_factor_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_element,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // evaluate
    test_basis.jacobians(point_in_reference_element, test_basis_grads_, param);
    ansatz_basis.jacobians(point_in_reference_element, ansatz_basis_grads_, param);
    const auto diffusion = local_diffusion_factor_->evaluate(point_in_reference_element, param);
    symmetric_ansatz_basis_grads_.resize(cols);
    for (size_t jj = 0; jj < cols; ++jj)
      for (size_t rr = 0; rr < r; ++rr)
        for (size_t cc = 0; cc < d; ++cc)
          symmetric_ansatz_basis_grads_[jj][rr][cc] =
              0.5 * (ansatz_basis_grads_[jj][rr][cc] + ansatz_basis_grads_[jj][cc][rr]);
    // compute elliptic evaluation
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        for (size_t rr = 0; rr < r; ++rr)
          for (size_t cc = 0; cc < d; ++cc)
            result[ii][jj] += symmetric_ansatz_basis_grads_[jj][rr][cc] * test_basis_grads_[ii][rr][cc] * diffusion;
  } // ... evaluate(...)

private:
  const std::unique_ptr<DiffusionFactorType> diffusion_factor_;
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> symmetric_ansatz_basis_grads_;
}; // class LocalSymmetrizedLaplaceIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_SYMMETRIZED_LAPLACE_HH
