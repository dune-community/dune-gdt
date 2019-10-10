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

#warning This header is deprecated, use and include <dune/gdt/local/integrands/laplace.hh> instead!

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
#  define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH

#  include <dune/xt/common/deprecated.hh>
#  include <dune/xt/common/memory.hh>
#  include <dune/xt/la/container/eye-matrix.hh>
#  include <dune/xt/functions/base/function-as-grid-function.hh>
#  include <dune/xt/functions/constant.hh>
#  include <dune/xt/functions/interfaces/grid-function.hh>

#  include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * This class is deprecated, use LocalLaplaceIntegrand instead (10.08.2019)!
 * Given an inducing scalar function lambda and an inducing matrix-valued function kappa, computes
 * `lambda(x) * {[kappa(x) \nabla phi(x)] * \nabla psi(x)}` for all combinations of phi in the ansatz basis and psi in
 * the test basis.
 * If phi and psi are vector-valued, \nabla phi is the jacobian matrix and we are actually computing
 * `lambda(x) * {[kappa(x) (\nabla phi(x))^T] : (\nabla psi(x))^T}`, where ':' denotes the matrix scalar product.
 */
template <class E, size_t r = 1, class F = double>
class /*DXT_DEPRECATED_MSG("Use LocalLaplaceIntegrand instead (10.08.2019)!")*/ LocalEllipticIntegrand
  : public LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, F, F, r, 1, F>;
  using ThisType = LocalEllipticIntegrand;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  LocalEllipticIntegrand(
      const F& diffusion_factor = F(1),
      const XT::Common::FieldMatrix<F, d, d>& diffusion_tensor = XT::LA::eye_matrix<FieldMatrix<F, d, d>>(d, d))
    : BaseType()
    , diffusion_factor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(
          new XT::Functions::ConstantFunction<d, 1, 1, F>(diffusion_factor)))
    , diffusion_tensor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, F>(
          new XT::Functions::ConstantFunction<d, d, d, F>(diffusion_tensor)))
    , local_diffusion_factor_(diffusion_factor_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const XT::Functions::FunctionInterface<d, 1, 1, F>& diffusion_factor,
                         const XT::Functions::FunctionInterface<d, d, d, F>& diffusion_tensor)
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , diffusion_factor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(diffusion_factor))
    , diffusion_tensor_(new XT::Functions::FunctionAsGridFunctionWrapper<E, d, d, F>(diffusion_tensor))
    , local_diffusion_factor_(diffusion_factor_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor)
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_(diffusion_factor_.access().local_function())
    , local_diffusion_tensor_(diffusion_tensor_.access().local_function())
  {}

  LocalEllipticIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_(diffusion_factor_.access().local_function())
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
#  ifndef NDEBUG
    if (!ele.geometry().affine())
      std::cerr << "Warning: integration order has to be increased for non-affine geometries!" << std::endl;
#  endif
    local_diffusion_factor_->bind(ele);
    local_diffusion_tensor_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_factor_->order(param) + local_diffusion_tensor_->order(param) + test_basis.order(param)
           + ansatz_basis.order(param);
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
  const XT::Common::ConstStorageProvider<DiffusionFactorType> diffusion_factor_;
  const XT::Common::ConstStorageProvider<DiffusionTensorType> diffusion_tensor_;
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class LocalEllipticIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_HH
