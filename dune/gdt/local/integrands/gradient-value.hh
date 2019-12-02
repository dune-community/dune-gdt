// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_GRADIENT_VALUE_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_GRADIENT_VALUE_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing (vector-valued) function u, computes `(\nabla phi(x) * u(x)) * psi(x)` for all combinations of phi
 * in the ansatz basis and psi in the test basis. If use_test_gradient is set to true, ansatz and test basis are
 * swapped.
 *
 * \sa local_binary_to_unary_element_integrand
 */
template <class E,
          size_t r = 1,
          size_t rC = 1,
          class TR = double,
          class F = double,
          class AR = TR,
          bool use_test_gradient = false>
class LocalElementGradientValueIntegrand : public LocalBinaryElementIntegrandInterface<E, r, rC, TR, F, r, rC, AR>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, rC, TR, F, r, rC, AR>;
  using ThisType = LocalElementGradientValueIntegrand;
  static_assert(rC == 1, "Not tested for rC > 1!");

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;
  using BasisValues =
      std::vector<typename std::conditional_t<use_test_gradient, LocalAnsatzBasisType, LocalTestBasisType>::RangeType>;
  using BasisJacobians = std::vector<
      typename std::conditional_t<use_test_gradient, LocalTestBasisType, LocalAnsatzBasisType>::DerivativeRangeType>;
  static const size_t vector_size = d;

  using VectorGridFunctionType = XT::Functions::GridFunction<E, vector_size, 1, F>;
  using VectorValues = typename VectorGridFunctionType::LocalFunctionType::RangeReturnType;

  LocalElementGradientValueIntegrand(const VectorGridFunctionType vector_in)
    : BaseType()
    , vector_(vector_in)
    , local_function_(vector_.local_function())
  {}

  LocalElementGradientValueIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , vector_(other.vector_)
    , local_function_(vector_.local_function())
  {}

  LocalElementGradientValueIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override final
  {
    local_function_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_function_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  using BaseType::evaluate;

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
    vector_values_ = local_function_->evaluate(point_in_reference_element, param);
    helper<>::evaluate(
        test_basis, ansatz_basis, point_in_reference_element, vector_values_, values_, jacobians_, result, param);
  } // ... evaluate(...)

private:
  template <bool use_test = use_test_gradient, bool anything = true>
  struct helper
  {
    static void evaluate(const LocalTestBasisType& test_basis,
                         const LocalAnsatzBasisType& ansatz_basis,
                         const DomainType& point_in_reference_element,
                         const VectorValues& vector_values,
                         BasisValues& values,
                         BasisJacobians& jacobians,
                         DynamicMatrix<F>& result,
                         const XT::Common::Parameter& param)
    {
      test_basis.evaluate(point_in_reference_element, values, param);
      ansatz_basis.jacobians(point_in_reference_element, jacobians, param);

      auto nabla_phi_times_u = values[0];
      const size_t rows = test_basis.size();
      const size_t cols = ansatz_basis.size();
      for (size_t jj = 0; jj < cols; ++jj) {
        jacobians[jj].mv(vector_values, nabla_phi_times_u);
        for (size_t ii = 0; ii < rows; ++ii)
          result[ii][jj] = nabla_phi_times_u * values[ii];
      }
    }
  };

  template <bool anything>
  struct helper<true, anything>
  {
    static void evaluate(const LocalTestBasisType& test_basis,
                         const LocalAnsatzBasisType& ansatz_basis,
                         const DomainType& point_in_reference_element,
                         const VectorValues& vector_values,
                         BasisValues& values,
                         BasisJacobians& jacobians,
                         DynamicMatrix<F>& result,
                         const XT::Common::Parameter& param)
    {
      ansatz_basis.evaluate(point_in_reference_element, values, param);
      test_basis.jacobians(point_in_reference_element, jacobians, param);

      auto nabla_psi_times_u = values[0];
      const size_t rows = test_basis.size();
      const size_t cols = ansatz_basis.size();
      for (size_t ii = 0; ii < rows; ++ii) {
        jacobians[ii].mv(vector_values, nabla_psi_times_u);
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] = nabla_psi_times_u * values[jj];
      }
    }
  };

  const VectorGridFunctionType vector_;
  std::unique_ptr<typename VectorGridFunctionType::LocalFunctionType> local_function_;
  mutable VectorValues vector_values_;
  mutable BasisValues values_;
  mutable BasisJacobians jacobians_;
}; // class LocalElementGradientValueIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_GRADIENT_VALUE_HH
