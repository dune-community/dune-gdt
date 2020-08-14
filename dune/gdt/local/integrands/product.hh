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
#include <dune/xt/grid/print.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/combined-functions.hh>
#include <dune/xt/functions/base/combined-grid-functions.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/print.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given an inducing function f (may be matrix-valued), computes `(f(x) * psi(x)) * phi(x)` for all combinations of phi
 * in the ansatz basis and psi in the test basis.
 *
 * \note Note that f can also be given as a scalar value or omitted.
 *
 * \note Applying f to the ansatz basis can be done by passing f^T (the transposed of f)
 *
 * \sa local_binary_to_unary_element_integrand
 */
template <class E, size_t r = 1, class TR = double, class F = double, class AR = TR>
class LocalElementProductIntegrand : public LocalBinaryElementIntegrandInterface<E, r, 1, TR, F, r, 1, AR>
{
  using ThisType = LocalElementProductIntegrand;
  using BaseType = LocalBinaryElementIntegrandInterface<E, r, 1, TR, F, r, 1, AR>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  LocalElementProductIntegrand(XT::Functions::GridFunction<E, r, r, F> weight = {1.},
                               const std::string& logging_prefix = "")
    : BaseType({},
               logging_prefix.empty() ? "ElementProductIntegrand" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , weight_(weight)
    , local_weight_(weight_.local_function())
  {
    LOG_(info) << "LocalElementProductIntegrand(this=" << this << ", weight=" << &weight << ")" << std::endl;
  }

  LocalElementProductIntegrand(const ThisType& other)
    : BaseType(other)
    , weight_(other.weight_)
    , local_weight_(weight_.local_function())
  {}

  LocalElementProductIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    LOG_(debug) << "copy_as_binary_element_integrand()" << std::endl;
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override final
  {
    LOG_(debug) << "post_bind(element=" << XT::Grid::print(ele) << ")" << std::endl;
    local_weight_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "order(element=" << XT::Grid::print(this->element()) << ", {test|ansatz}_basis.size()={"
                << test_basis.size(param) << "|" << ansatz_basis.size(param) << "}, param=" << param
                << ")\n     local_weight.order() = " << local_weight_->order(param)
                << "\n     test_basis.order() = " << test_basis.order(param)
                << "\n     ansatz_basis.order() = " << ansatz_basis.order(param) << "\n     returning "
                << local_weight_->order(param) + test_basis.order(param) + ansatz_basis.order(param) << std::endl;
    return local_weight_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_element,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "evaluate({test|ansatz}_basis.size()={" << test_basis.size(param) << "|" << ansatz_basis.size(param)
                << "}, point_in_{reference_element|physical_space} = {" << print(point_in_reference_element) << "|"
                << print(this->element().geometry().global(point_in_reference_element)) << "}, param=" << param << ")"
                << std::endl;
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    // evaluate
    test_basis.evaluate(point_in_reference_element, test_basis_values_, param);
    ansatz_basis.evaluate(point_in_reference_element, ansatz_basis_values_, param);
    const auto weight = local_weight_->evaluate(point_in_reference_element, param);
    LOG_(debug) << "  test_basis_values_ = " << test_basis_values_
                << "\n  ansatz_basis_values_ = " << ansatz_basis_values_ << "\n  weight = " << weight << std::endl;
    // compute product
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] = (weight * test_basis_values_[ii]) * ansatz_basis_values_[jj];
    LOG_(debug) << "  result = " << print(result, {{"oneline", "true"}}) << std::endl;
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, r, r, F> weight_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, r, r, F>::LocalFunctionType> local_weight_;
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
class LocalCouplingIntersectionProductIntegrand
  : public LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>
{
  using ThisType = LocalCouplingIntersectionProductIntegrand;
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GridFunctionType = XT::Functions::GridFunctionInterface<E, r, r, F>;

  LocalCouplingIntersectionProductIntegrand(XT::Functions::GridFunction<E, r, r, F> weight = {1.},
                                            const std::string& logging_prefix = "")
    : BaseType({},
               logging_prefix.empty() ? "gdt" : "gdt.localcouplingintersectionproductintegrand",
               logging_prefix.empty() ? "LocalCouplingIntersectionProductIntegrand" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , weight_(weight)
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  LocalCouplingIntersectionProductIntegrand(const ThisType& other)
    : BaseType(other)
    , weight_(other.weight_)
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  LocalCouplingIntersectionProductIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersct) override final
  {
    auto inside_element = intersct.inside();
    local_weight_in_->bind(inside_element);
    if (intersct.neighbor()) {
      local_weight_out_->bind(intersct.outside());
    } else
      local_weight_out_->bind(intersct.inside());
  } // ... post_bind(...)

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(local_weight_in_->order(param), local_weight_out_->order(param))
           + std::max(test_basis_inside.order(param), test_basis_outside.order(param))
           + std::max(ansatz_basis_inside.order(param), ansatz_basis_outside.order(param));
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
    const auto weight_in = local_weight_in_->evaluate(point_in_inside_reference_element, param);
    const auto point_in_outside_reference_element =
        this->intersection().geometryInOutside().global(point_in_reference_intersection);
    test_basis_outside.evaluate(point_in_outside_reference_element, test_basis_out_values_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    const auto weight_out = local_weight_out_->evaluate(point_in_outside_reference_element, param);
    // compute integrand
    const auto average_function_value = (weight_in + weight_out) * 0.5;
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
  XT::Functions::GridFunction<E, r, r, F> weight_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, r, r, F>::LocalFunctionType> local_weight_in_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, r, r, F>::LocalFunctionType> local_weight_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; // class LocalCouplingIntersectionProductIntegrand


/**
 * Given an inducing function f, computes `f * phi * psi` for all combinations of phi and psi in the bases.
 *
 * \note Note that f can also be given as a scalar value or omitted.
 */
template <class I, size_t r = 1, class TR = double, class F = double, class AR = TR>
class LocalIntersectionProductIntegrand : public LocalBinaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>
{
  using ThisType = LocalIntersectionProductIntegrand;
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, r, 1, TR, F, r, 1, AR>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GridFunctionType = XT::Functions::GridFunctionInterface<E, r, r, F>;

  LocalIntersectionProductIntegrand(XT::Functions::GridFunction<E, r, r, F> weight = {1.},
                                    const bool use_inside_bases = true,
                                    const std::string& logging_prefix = "")
    : BaseType({},
               logging_prefix.empty() ? "gdt" : "gdt.localintersectionproductintegrand",
               logging_prefix.empty() ? "LocalIntersectionProductIntegrand" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , weight_(weight)
    , local_weight_(weight_.local_function())
    , inside_(use_inside_bases)
  {}

  LocalIntersectionProductIntegrand(const ThisType& other)
    : BaseType(other)
    , weight_(other.weight_)
    , local_weight_(weight_.local_function())
    , inside_(other.inside_)
  {}

  LocalIntersectionProductIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersct) override final
  {
    if (inside_)
      local_weight_->bind(intersct.inside());
    else {
      DUNE_THROW_IF(!intersct.neighbor(), Exceptions::integrand_error, "");
      local_weight_->bind(intersct.outside());
    }
  } // ... post_bind(...)

public:
  bool inside() const override final
  {
    return inside_;
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_weight_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare sotrage
    this->ensure_size_and_clear_results(test_basis, ansatz_basis, result, param);
    // evaluate
    typename XT::Functions::GridFunction<E, r, r, F>::LocalFunctionType::RangeReturnType weight;
    if (inside_) {
      const auto point_in_inside_reference_element =
          this->intersection().geometryInInside().global(point_in_reference_intersection);
      test_basis.evaluate(point_in_inside_reference_element, test_basis_values_, param);
      ansatz_basis.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
      weight = local_weight_->evaluate(point_in_inside_reference_element, param);
    } else {
      const auto point_in_outside_reference_element =
          this->intersection().geometryInOutside().global(point_in_reference_intersection);
      test_basis.evaluate(point_in_outside_reference_element, test_basis_values_, param);
      ansatz_basis.evaluate(point_in_outside_reference_element, ansatz_basis_values_, param);
      weight = local_weight_->evaluate(point_in_outside_reference_element, param);
    }
    // compute integrand
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] = (weight * ansatz_basis_values_[jj]) * test_basis_values_[ii];
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, r, r, F> weight_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, r, r, F>::LocalFunctionType> local_weight_;
  const bool inside_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // class LocalIntersectionProductIntegrand


template <class E_or_I, size_t r = 1, class TR = double, class F = double, class AR = TR>
class LocalProductIntegrand
  : public std::conditional_t<XT::Grid::is_entity<E_or_I>::value,
                              LocalElementProductIntegrand<E_or_I, r, TR, F, AR>,
                              LocalIntersectionProductIntegrand<E_or_I, r, TR, F, AR>>
{
  using BaseType = std::conditional_t<XT::Grid::is_entity<E_or_I>::value,
                                      LocalElementProductIntegrand<E_or_I, r, TR, F, AR>,
                                      LocalIntersectionProductIntegrand<E_or_I, r, TR, F, AR>>;

public:
  template <class... Args>
  explicit LocalProductIntegrand(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
}; // class LocalProductIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_PRODUCTS_HH
