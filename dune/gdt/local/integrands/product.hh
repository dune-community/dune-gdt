// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2014, 2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/combined-functions.hh>
#include <dune/xt/functions/base/combined-grid-functions.hh>
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
template <class E, size_t r = 1, class TR = double, class F = double, class AR = TR>
class LocalElementProductIntegrand : public LocalBinaryElementIntegrandInterface<E, r, 1, TR, F, r, 1, AR>
{
  using ThisType = LocalElementProductIntegrand<E, r, TR, F, AR>;
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, TR, F, r, 1, AR>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GridFunctionType = XT::Functions::GridFunctionInterface<E, r, r, F>;

  LocalElementProductIntegrand(const F& inducing_value = F(1))
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductFunction<XT::Functions::ConstantFunction<d, 1, 1, F>,
                                             XT::Functions::ConstantFunction<d, r, r, F>>(
              new XT::Functions::ConstantFunction<d, 1, 1, F>(inducing_value),
              new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r)))))
    , local_function_(inducing_function_.access().local_function())
  {}

  LocalElementProductIntegrand(const XT::Functions::FunctionInterface<d, 1, 1, F>& inducing_function)
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductFunction<XT::Functions::FunctionInterface<d, 1, 1, F>,
                                             XT::Functions::ConstantFunction<d, r, r, F>>(
              inducing_function,
              new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r)))))
    , local_function_(inducing_function_.access().local_function())
  {}

  LocalElementProductIntegrand(const XT::Functions::GridFunctionInterface<E, 1, 1, F>& inducing_function)
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductGridFunction<XT::Functions::GridFunctionInterface<E, 1, 1, F>,
                                                 XT::Functions::GridFunctionInterface<E, r, r, F>>(
              inducing_function,
              new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
                  new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r))))))
    , local_function_(inducing_function_.access().local_function())
  {}

  template <class E_, typename = std::enable_if_t<std::is_same<E_, E>::value && r != 1, void>>
  LocalElementProductIntegrand(const XT::Functions::GridFunctionInterface<E_, r, r, F>& inducing_function)
    : BaseType()
    , inducing_function_(inducing_function)
    , local_function_(inducing_function_.access().local_function())
  {}

  LocalElementProductIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , inducing_function_(other.inducing_function_)
    , local_function_(inducing_function_.access().local_function())
  {}

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
        result[ii][jj] = (function_value * test_basis_values_[ii]) * ansatz_basis_values_[jj];
  } // ... evaluate(...)

private:
  const XT::Common::ConstStorageProvider<GridFunctionType> inducing_function_;
  std::unique_ptr<typename GridFunctionType::LocalFunctionType> local_function_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // class LocalElementProductIntegrand


/**
 * Given an inducing function f, computes `<f> * phi * psi` for all combinations of phi and psi in the bases, where
 * `<f>` denotes the average of f evaluated on the inside and evaluated on the outside.
 *
 * \note Note that f can also be given as a scalar value or omitted.
 */
template <class I, size_t r = 1, class TR = double, class F = double, class AR = TR>
class LocalIntersectionProductIntegrand : public LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>
{
  using ThisType = LocalIntersectionProductIntegrand<I, r, TR, F, AR>;
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GridFunctionType = XT::Functions::GridFunctionInterface<E, r, r, F>;

  LocalIntersectionProductIntegrand(const F& inducing_value = F(1))
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductFunction<XT::Functions::ConstantFunction<d, 1, 1, F>,
                                             XT::Functions::ConstantFunction<d, r, r, F>>(
              new XT::Functions::ConstantFunction<d, 1, 1, F>(inducing_value),
              new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r)))))
    , local_function_in_(inducing_function_.access().local_function())
    , local_function_out_(inducing_function_.access().local_function())
  {}

  LocalIntersectionProductIntegrand(const XT::Functions::FunctionInterface<d, 1, 1, F>& inducing_function)
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductFunction<XT::Functions::FunctionInterface<d, 1, 1, F>,
                                             XT::Functions::ConstantFunction<d, r, r, F>>(
              inducing_function,
              new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r)))))
    , local_function_in_(inducing_function_.access().local_function())
    , local_function_out_(inducing_function_.access().local_function())
  {}

  LocalIntersectionProductIntegrand(const XT::Functions::GridFunctionInterface<E, 1, 1, F>& inducing_function)
    : BaseType()
    , inducing_function_(new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
          new XT::Functions::ProductGridFunction<XT::Functions::GridFunctionInterface<E, 1, 1, F>,
                                                 XT::Functions::GridFunctionInterface<E, r, r, F>>(
              inducing_function,
              new XT::Functions::FunctionAsGridFunctionWrapper<E, r, r, F>(
                  new XT::Functions::ConstantFunction<d, r, r, F>(XT::LA::eye_matrix<FieldMatrix<F, r, r>>(r))))))
    , local_function_in_(inducing_function_.access().local_function())
    , local_function_out_(inducing_function_.access().local_function())
  {}

  template <class E_, typename = std::enable_if_t<std::is_same<E_, E>::value && r != 1, void>>
  LocalIntersectionProductIntegrand(const XT::Functions::GridFunctionInterface<E_, r, r, F>& inducing_function)
    : BaseType()
    , inducing_function_(inducing_function)
    , local_function_in_(inducing_function_.access().local_function())
    , local_function_out_(inducing_function_.access().local_function())
  {}

  LocalIntersectionProductIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , inducing_function_(other.inducing_function_)
    , local_function_in_(inducing_function_.access().local_function())
    , local_function_out_(inducing_function_.access().local_function())
  {}

  LocalIntersectionProductIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersct) override final
  {
    auto inside_element = intersct.inside();
    local_function_in_->bind(inside_element);
    if (intersct.neighbor()) {
      local_function_out_->bind(intersct.outside());
    } else
      local_function_out_->bind(intersct.inside());
  } // ... post_bind(...)

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(local_function_in_->order(param), local_function_out_->order(param))
           + std::max(test_basis_inside.order(param), test_basis_outside.order(param))
           + std::max(ansatz_basis_outside.order(param), ansatz_basis_outside.order(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis_inside,
                const LocalAnsatzBasisType& ansatz_basis_inside,
                const LocalTestBasisType& test_basis_outside,
                const LocalAnsatzBasisType& ansatz_basis_outside,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result_in_in,
                DynamicMatrix<F>& result_in_out,
                DynamicMatrix<F>& result_out_in,
                DynamicMatrix<F>& result_out_out,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare sotrage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size = [](auto& m, const auto& rws, const auto& cls) {
      if (m.rows() < rws || m.cols() < cls)
        m.resize(rws, cls);
    };
    ensure_size(result_in_in, rows_in, cols_in);
    ensure_size(result_in_out, rows_in, cols_out);
    ensure_size(result_out_in, rows_out, cols_in);
    ensure_size(result_out_out, rows_out, cols_out);
    // evaluate
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_in_values_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_in_values_, param);
    const auto function_value_in = local_function_in_->evaluate(point_in_inside_reference_element, param);
    const auto point_in_outside_reference_element =
        this->intersection().geometryInOutside().global(point_in_reference_intersection);
    test_basis_outside.evaluate(point_in_outside_reference_element, test_basis_out_values_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    const auto function_value_out = local_function_out_->evaluate(point_in_outside_reference_element, param);
    // compute integrand
    const auto average_function_value = (function_value_in + function_value_out) * 0.5;
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_in_in[ii][jj] = (average_function_value * ansatz_basis_in_values_[jj]) * test_basis_in_values_[ii];
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_in_out[ii][jj] = (average_function_value * ansatz_basis_out_values_[jj]) * test_basis_in_values_[ii];
    }
    for (size_t ii = 0; ii < rows_out; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_out_in[ii][jj] = (average_function_value * ansatz_basis_in_values_[jj]) * test_basis_out_values_[ii];
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_out_out[ii][jj] = (average_function_value * ansatz_basis_out_values_[jj]) * test_basis_out_values_[ii];
    }
  } // ... evaluate(...)

private:
  const XT::Common::ConstStorageProvider<GridFunctionType> inducing_function_;
  std::unique_ptr<typename GridFunctionType::LocalFunctionType> local_function_in_;
  std::unique_ptr<typename GridFunctionType::LocalFunctionType> local_function_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; // class LocalIntersectionProductIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
