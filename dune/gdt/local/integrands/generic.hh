// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_GENERIC_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_GENERIC_HH

#include <functional>

#include <dune/gdt/exceptions.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r = 1, size_t rC = 1, class R = double, class F = double>
class GenericLocalUnaryElementIntegrand : public LocalUnaryElementIntegrandInterface<E, r, rC, R, F>
{
  using ThisType = GenericLocalUnaryElementIntegrand;
  using BaseType = LocalUnaryElementIntegrandInterface<E, r, rC, R, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType =
      std::function<int(const LocalTestBasisType& /*basis*/, const XT::Common::Parameter& /*param*/)>;
  using GenericEvaluateFunctionType = std::function<void(const LocalTestBasisType& /*basis*/,
                                                         const DomainType& /*point_in_reference_element*/,
                                                         DynamicVector<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;
  using GenericPostBindFunctionType = std::function<void(const E& /*ele*/)>;

  GenericLocalUnaryElementIntegrand(
      GenericOrderFunctionType order_function,
      GenericEvaluateFunctionType evaluate_function,
      GenericPostBindFunctionType post_bind_function = [](const E&) {},
      const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
    , post_bind_(post_bind_function)
  {}

  GenericLocalUnaryElementIntegrand(const ThisType& other)
    : BaseType(other)
    , order_(other.order_)
    , evaluate_(other.evaluate_)
    , post_bind_(other.post_bind_)
  {}

  std::unique_ptr<BaseType> copy_as_unary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const E& ele) override final
  {
    post_bind_(ele);
  }

public:
  int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const override final
  {
    return order_(basis, this->parse_parameter(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& basis,
                const DomainType& point_in_reference_element,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    // evaluate
    evaluate_(basis, point_in_reference_element, result, this->parse_parameter(param));
    // check
    DUNE_THROW_IF(result.size() < size,
                  Exceptions::integrand_error,
                  "basis.size(param) = " << size << "\n   result.size() = " << result.size());
  } // ... evaluate(...)

private:
  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericPostBindFunctionType post_bind_;
}; // class GenericLocalUnaryElementIntegrand


template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class GenericLocalBinaryElementIntegrand
  : public LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using ThisType = GenericLocalBinaryElementIntegrand;
  using BaseType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType = std::function<int(const LocalTestBasisType& /*test_basis*/,
                                                     const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                     const XT::Common::Parameter& /*param*/)>;
  using GenericEvaluateFunctionType = std::function<void(const LocalTestBasisType& /*test_basis*/,
                                                         const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                         const DomainType& /*point_in_reference_element*/,
                                                         DynamicMatrix<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;
  using GenericPostBindFunctionType = std::function<void(const E& /*ele*/)>;

  GenericLocalBinaryElementIntegrand(
      GenericOrderFunctionType order_function,
      GenericEvaluateFunctionType evaluate_function,
      GenericPostBindFunctionType post_bind_function = [](const E&) {},
      const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
    , post_bind_(post_bind_function)
  {}

  GenericLocalBinaryElementIntegrand(const ThisType& other)
    : BaseType(other)
    , order_(other.order_)
    , evaluate_(other.evaluate_)
    , post_bind_(other.post_bind_)
  {}

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const E& ele) override final
  {
    post_bind_(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return order_(test_basis, ansatz_basis, this->parse_parameter(param));
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
    evaluate_(test_basis, ansatz_basis, point_in_reference_element, result, this->parse_parameter(param));
    // check
    DUNE_THROW_IF(result.rows() < rows || result.cols() < cols,
                  Exceptions::integrand_error,
                  "test_basis.size(param) = " << rows << "\n   result.rows() = " << result.rows()
                                              << "ansatz_basis.size(param) = " << cols
                                              << "\n   result.cols() = " << result.cols());
  } // ... evaluate(...)

private:
  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericPostBindFunctionType post_bind_;
}; // class GenericLocalBinaryElementIntegrand


template <class I, size_t r = 1, size_t rC = 1, class RF = double, class F = double>
class GenericLocalUnaryIntersectionIntegrand : public LocalUnaryIntersectionIntegrandInterface<I, r, rC, RF, F>
{
  using ThisType = GenericLocalUnaryIntersectionIntegrand;
  using BaseType = LocalUnaryIntersectionIntegrandInterface<I, r, rC, RF, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType =
      std::function<int(const LocalTestBasisType& /*basis*/, const XT::Common::Parameter& /*param*/)>;
  using GenericEvaluateFunctionType = std::function<void(const LocalTestBasisType& /*basis*/,
                                                         const DomainType& /*point_in_reference_intersection*/,
                                                         DynamicVector<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;
  using GenericPostBindFunctionType = std::function<void(const I& /*intersection*/)>;

  GenericLocalUnaryIntersectionIntegrand(
      GenericOrderFunctionType order_function,
      GenericEvaluateFunctionType evaluate_function,
      GenericPostBindFunctionType post_bind_function = [](const I&) {},
      bool use_inside_bases = true,
      const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
    , post_bind_(post_bind_function)
    , inside_(use_inside_bases)
  {}

  GenericLocalUnaryIntersectionIntegrand(const ThisType& other) = default;

  std::unique_ptr<BaseType> copy_as_unary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const I& intersect) override final
  {
    post_bind_(intersect);
  }

public:
  bool inside() const override final
  {
    return inside_;
  }

  int order(const LocalTestBasisType& basis, const XT::Common::Parameter& param = {}) const override final
  {
    return order_(basis, this->parse_parameter(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& basis,
                const DomainType& point_in_reference_intersection,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    // evaluate
    evaluate_(basis, point_in_reference_intersection, result, this->parse_parameter(param));
    // check
    DUNE_THROW_IF(result.size() < size,
                  Exceptions::integrand_error,
                  "basis.size(param) = " << size << "\n   result.size() = " << result.size());
  } // ... evaluate(...)

private:
  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericPostBindFunctionType post_bind_;
  const bool inside_;
}; // class GenericLocalUnaryIntersectionIntegrand


template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class GenericLocalBinaryIntersectionIntegrand
  : public LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>
{
  using ThisType = GenericLocalBinaryIntersectionIntegrand;
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType = std::function<int(const LocalTestBasisType& /*test_basis*/,
                                                     const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                     const XT::Common::Parameter& /*param*/)>;
  using GenericEvaluateFunctionType = std::function<void(const LocalTestBasisType& /*test_basis*/,
                                                         const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                         const DomainType& /*point_in_reference_intersection*/,
                                                         DynamicMatrix<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;
  using GenericPostBindFunctionType = std::function<void(const I& /*intersection*/)>;

  GenericLocalBinaryIntersectionIntegrand(
      GenericOrderFunctionType order_function,
      GenericEvaluateFunctionType evaluate_function,
      GenericPostBindFunctionType post_bind_function = [](const I&) {},
      bool use_inside_bases = true,
      const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
    , post_bind_(post_bind_function)
    , inside_(use_inside_bases)
  {}

  GenericLocalBinaryIntersectionIntegrand(const ThisType& other) = default;

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  void post_bind(const I& intersect) override final
  {
    post_bind_(intersect);
  }

  bool inside() const override final
  {
    return inside_;
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return order_(test_basis, ansatz_basis, this->parse_parameter(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    // evaluate
    evaluate_(test_basis, ansatz_basis, point_in_reference_intersection, result, this->parse_parameter(param));
    // check
    DUNE_THROW_IF(result.rows() < rows || result.cols() < cols,
                  Exceptions::integrand_error,
                  "test_basis.size(param) = " << rows << "\n   result.rows() = " << result.rows()
                                              << "ansatz_basis.size(param) = " << cols
                                              << "\n   result.cols() = " << result.cols());
  } // ... evaluate(...)

private:
  const GenericOrderFunctionType order_;
  const GenericEvaluateFunctionType evaluate_;
  const GenericPostBindFunctionType post_bind_;
  const bool inside_;
}; // class GenericLocalBinaryIntersectionIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_GENERIC_HH
