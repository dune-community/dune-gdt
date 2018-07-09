// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing scalar function lambda and an inducing matrix-valued function kappa, computes
 * `lambda(x) * {[kappa(x) \nabla phi(x)] * \nabla psi(x)}` for all combinations of phi and psi in the bases.
   */
template <class E, size_t r = 1, class F = double>
class LocalEllipticIntegrand : public LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using ThisType = LocalEllipticIntegrand<E, r, F>;

public:
  using BaseType::d;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;
  using typename BaseType::DomainType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  LocalEllipticIntegrand(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor)
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {
  }

  LocalEllipticIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {
  }

  LocalEllipticIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele)
  {
    local_diffusion_factor_->bind(ele);
    local_diffusion_tensor_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_factor_->order(param) + local_diffusion_tensor_->order(param)
           + std::max(test_basis.order(param) - 1, 0) + std::max(ansatz_basis.order(param) - 1, 0);
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
    const auto diffusion = local_diffusion_tensor_->evaluate(point_in_reference_element, param)
                           * local_diffusion_factor_->evaluate(point_in_reference_element, param);
    // compute elliptic evaluation
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        for (size_t rr = 0; rr < r; ++rr)
          result[ii][jj] += (diffusion * ansatz_basis_grads_[jj][rr]) * test_basis_grads_[ii][rr];
  } // ... evaluate(...)

private:
  const DiffusionFactorType& diffusion_factor_; // These are just required ...
  const DiffusionTensorType& diffusion_tensor_; //                         ... for the copy ctor atm.
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class LocalEllipticIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
