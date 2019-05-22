// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given a weight function kappa, computes `{[kappa(x) \nabla phi(x)] * \nabla psi(x)}` for all combinations of phi and
 * psi in the bases.
 */
template <class E, size_t r = 1, class F = double>
class LocalLaplaceIntegrand : public LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using ThisType = LocalLaplaceIntegrand;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  explicit LocalLaplaceIntegrand(
      XT::Functions::GridFunction<E, d, d, F> diffusion = XT::LA::eye_matrix<FieldMatrix<F, d, d>>(d, d))
    : BaseType()
    , weight_(diffusion)
    , local_weight_(weight_.local_function())
  {}

  LocalLaplaceIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , weight_(other.weight_)
    , local_weight_(weight_.local_function())
  {}

  LocalLaplaceIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override
  {
#ifndef NDEBUG
    if (!ele.geometry().affine())
      std::cerr << "Warning: integration order has to be increased for non-affine geometries!" << std::endl;
#endif
    local_weight_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_weight_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
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
    const auto weight = local_weight_->evaluate(point_in_reference_element, param);
    // compute elliptic evaluation
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        for (size_t rr = 0; rr < r; ++rr)
          result[ii][jj] += (weight * ansatz_basis_grads_[jj][rr]) * test_basis_grads_[ii][rr];
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, d, d, F> weight_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, d, F>::LocalFunctionType> local_weight_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class LocalLaplaceIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_HH
