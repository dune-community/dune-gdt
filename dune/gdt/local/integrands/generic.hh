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
  using ThisType = GenericLocalUnaryElementIntegrand<E, r, rC, R, F>;
  using BaseType = LocalUnaryElementIntegrandInterface<E, r, rC, R, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalBasisType;

  using GenericOrderFunctionType =
      std::function<int(const LocalBasisType& /*basis*/, const XT::Common::Parameter& /*param*/)>;
  using GenericEvalauteFunctionType = std::function<void(const LocalBasisType& /*basis*/,
                                                         const DomainType& /*point_in_reference_element*/,
                                                         DynamicVector<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;

  GenericLocalUnaryElementIntegrand(GenericOrderFunctionType order_function,
                                    GenericEvalauteFunctionType evaluate_function,
                                    const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
  {}

  GenericLocalUnaryElementIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , order_(other.order_)
    , evaluate_(other.evaluate_)
  {}

  std::unique_ptr<BaseType> copy_as_unary_element_integrand() const
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const
  {
    return order_(basis, this->parse_parameter(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalBasisType& basis,
                const DomainType& point_in_reference_element,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const
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
  const GenericEvalauteFunctionType evaluate_;
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
  using ThisType = GenericLocalBinaryElementIntegrand<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using BaseType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType = std::function<int(const LocalTestBasisType& /*test_basis*/,
                                                     const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                     const XT::Common::Parameter& /*param*/)>;
  using GenericEvalauteFunctionType = std::function<void(const LocalTestBasisType& /*test_basis*/,
                                                         const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                         const DomainType& /*point_in_reference_element*/,
                                                         DynamicMatrix<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;

  GenericLocalBinaryElementIntegrand(GenericOrderFunctionType order_function,
                                     GenericEvalauteFunctionType evaluate_function,
                                     const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
  {}

  GenericLocalBinaryElementIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , order_(other.order_)
    , evaluate_(other.evaluate_)
  {}

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const
  {
    return order_(test_basis, ansatz_basis, this->parse_parameter(param));
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_element,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const
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
  const GenericEvalauteFunctionType evaluate_;
}; // class GenericLocalBinaryElementIntegrand


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
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using ThisType = GenericLocalBinaryIntersectionIntegrand<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericOrderFunctionType = std::function<int(const LocalTestBasisType& /*test_basis*/,
                                                     const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                     const XT::Common::Parameter& /*param*/)>;
  using GenericEvalauteFunctionType = std::function<void(const LocalTestBasisType& /*test_basis*/,
                                                         const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                         const DomainType& /*point_in_reference_intersection*/,
                                                         DynamicMatrix<F>& /*result*/,
                                                         const XT::Common::Parameter& /*param*/)>;

  GenericLocalBinaryIntersectionIntegrand(GenericOrderFunctionType order_function,
                                          GenericEvalauteFunctionType evaluate_function,
                                          bool use_inside_bases = true,
                                          const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , order_(order_function)
    , evaluate_(evaluate_function)
    , inside_(use_inside_bases)
  {}

  GenericLocalBinaryIntersectionIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , order_(other.order_)
    , evaluate_(other.evaluate_)
    , inside_(other.inside_)
  {}

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
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
  { // prepare storage
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
  const GenericEvalauteFunctionType evaluate_;
  const bool inside_;
}; // class GenericLocalBinaryIntersectionIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_GENERIC_HH
