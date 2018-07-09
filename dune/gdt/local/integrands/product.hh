// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing function f, computes `f(x) * phi(x) * psi(x)` for all combinations of phi and psi in the bases.
 *
 * \note Note that f can also be given as a scalar value or omitted.
 *
 * \sa local_binary_to_unary_element_integrand
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalElementProductIntegrand : public LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = LocalElementProductIntegrand<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;
  using BaseType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using typename BaseType::DomainType;
  using BaseType::d;
  using typename BaseType::ElementType;
  using typename BaseType::LocalTestBasisType;
  using typename BaseType::LocalAnsatzBasisType;

  using GridFunctionType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;

  LocalElementProductIntegrand(const F& inducing_value = F(1))
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(
          new XT::Functions::ConstantFunction<d, 1, 1, F>(inducing_value)))
    , local_function_(inducing_function_.access().local_function())
    , test_basis_values_()
    , ansatz_basis_values_()
  {
  }

  LocalElementProductIntegrand(const XT::Functions::FunctionInterface<d, 1, 1, F>& inducing_function)
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, F>(inducing_function))
    , local_function_(inducing_function_.access().local_function())
    , test_basis_values_()
    , ansatz_basis_values_()
  {
  }

  LocalElementProductIntegrand(const GridFunctionType& inducing_function)
    : BaseType(inducing_function.parameter_type())
    , inducing_function_(inducing_function)
    , local_function_(inducing_function_.access().local_function())
    , test_basis_values_()
    , ansatz_basis_values_()
  {
  }

  LocalElementProductIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , inducing_function_(other.inducing_function_)
    , local_function_(inducing_function_.access().local_function())
    , test_basis_values_()
    , ansatz_basis_values_()
  {
  }

  LocalElementProductIntegrand(ThisType&& source) = default;

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
    // evaluate
    test_basis.evaluate(point_in_reference_element, test_basis_values_, param);
    ansatz_basis.evaluate(point_in_reference_element, ansatz_basis_values_, param);
    const auto function_value = local_function_->evaluate(point_in_reference_element, param);
    // compute product
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] = function_value * test_basis_values_[ii] * ansatz_basis_values_[jj];
  } // ... evaluate(...)

private:
  const XT::Common::ConstStorageProvider<GridFunctionType> inducing_function_;
  std::unique_ptr<typename GridFunctionType::LocalFunctionType> local_function_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // class LocalElementProductIntegrand


/// \todo add LocalIntersectionProductIntegrand
///

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
