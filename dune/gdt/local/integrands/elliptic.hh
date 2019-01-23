// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   Kirsten Weber   (2013)
//   René Fritze     (2014, 2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016 - 2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing scalar function lambda and an inducing matrix-valued function kappa, computes
 * `lambda(x) * {[kappa(x) \nabla phi(x)] * \nabla psi(x)}` for all combinations of phi and psi in the bases.
 * If symmetrize is set to true, uses the symmetric part of the ansatz gradient, i.e., in the formula above
 * \nabla \phi(x) is replaced by 0.5 (\nabla \phi(x) + \nabla \phi(x)^T)
 */
template <class E, size_t r = 1, class F = double, bool symmetrize = false>
class LocalEllipticIntegrand : public LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using ThisType = LocalEllipticIntegrand;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;
  static_assert(!symmetrize || d == r, "To use the symmetric part the gradient has to be a square matrix!");

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  LocalEllipticIntegrand(
      const F& inducing_function = F(1),
      const XT::Common::FieldMatrix<F, d, d>& diffusion_tensor = XT::LA::eye_matrix<FieldMatrix<F, d, d>>(d, d))
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(
          new XT::Functions::ConstantFunction<d, 1, 1, F>(inducing_function)))
    , diffusion_tensor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, F>(
          new XT::Functions::ConstantFunction<d, d, d, F>(diffusion_tensor)))
    , local_inducing_function_(inducing_function_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const XT::Functions::FunctionInterface<d, 1, 1, F>& inducing_function,
                         const XT::Functions::FunctionInterface<d, d, d, F>& diffusion_tensor)
    : BaseType(inducing_function.parameter_type() + diffusion_tensor.parameter_type())
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(inducing_function))
    , diffusion_tensor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, F>(diffusion_tensor))
    , local_inducing_function_(inducing_function_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const DiffusionFactorType& inducing_function, const DiffusionTensorType& diffusion_tensor)
    : BaseType(inducing_function.parameter_type() + diffusion_tensor.parameter_type())
    , inducing_function_(inducing_function)
    , diffusion_tensor_(diffusion_tensor)
    , local_inducing_function_(inducing_function_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , inducing_function_(other.inducing_function_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_inducing_function_(inducing_function_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override
  {
    local_inducing_function_->bind(ele);
    local_diffusion_tensor_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_inducing_function_->order(param) + local_diffusion_tensor_->order(param)
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
                           * local_inducing_function_->evaluate(point_in_reference_element, param);
    // compute elliptic evaluation
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        for (size_t rr = 0; rr < r; ++rr)
          result[ii][jj] += (diffusion * ansatz_basis_grads_[jj][rr]) * test_basis_grads_[ii][rr];
    if (symmetrize) {
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          for (size_t rr = 0; rr < r; ++rr)
            for (size_t dd = 0; dd < d; ++dd)
              result[ii][jj] += (diffusion * ansatz_basis_grads_[jj][rr][dd]) * test_basis_grads_[ii][dd][rr];
      result *= 0.5;
    }
  } // ... evaluate(...)

private:
  const XT::Common::ConstStorageProvider<DiffusionFactorType> inducing_function_;
  const XT::Common::ConstStorageProvider<DiffusionTensorType> diffusion_tensor_;
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_inducing_function_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class LocalEllipticIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
