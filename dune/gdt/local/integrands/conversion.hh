// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH

#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/functions/generic/element-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given a function f and a binary element integrand bi(\cdot, \cdot), models the unary element integrand bi(f, \cdot).
 *
 * Most likey you do not want to use this class directly, but LocalUnaryElementIntegrandInterface::with_ansatz.
 *
 * See also LocalUnaryElementIntegrandInterface for a description of the template arguments.
 *
 * \sa LocalUnaryElementIntegrandInterface
 * \sa LocalBinaryElementIntegrandInterface
 * \sa XT::Functions::GridFunction
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryToUnaryElementIntegrand : public LocalUnaryElementIntegrandInterface<E, t_r, t_rC, TF, F>
{
  using ThisType = LocalBinaryToUnaryElementIntegrand;
  using BaseType = LocalUnaryElementIntegrandInterface<E, t_r, t_rC, TF, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalBasisType;

  using LocalBinaryElementIntegrandType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;

  LocalBinaryToUnaryElementIntegrand(const LocalBinaryElementIntegrandType& local_binary_integrand,
                                     XT::Functions::GridFunction<E, a_r, a_rC, AF> inducing_function_as_ansatz_basis,
                                     const std::string& logging_prefix = "")
    : BaseType({},
               logging_prefix.empty() ? "gdt" : "gdt.localbinarytounaryelementintegrand",
               logging_prefix.empty() ? "LocalBinaryToUnaryElementIntegrand" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , inducing_function_as_ansatz_basis_(inducing_function_as_ansatz_basis)
    , local_function_(inducing_function_as_ansatz_basis_.local_function())
    , local_binary_integrand_(local_binary_integrand.copy_as_binary_element_integrand())
  {
    LOG_(debug) << this->logging_id << "(local_binary_integrand=" << &local_binary_integrand
                << ", inducing_function_as_ansatz_basis=" << &inducing_function_as_ansatz_basis << ")" << std::endl;
  }

  LocalBinaryToUnaryElementIntegrand(const ThisType& other)
    : BaseType(other)
    , inducing_function_as_ansatz_basis_(other.inducing_function_as_ansatz_basis_)
    , local_function_(inducing_function_as_ansatz_basis_.local_function())
    , local_binary_integrand_(other.local_binary_integrand_->copy_as_binary_element_integrand())
  {}

  std::unique_ptr<BaseType> copy_as_unary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& element) override final
  {
    local_function_->bind(element);
    local_binary_integrand_->bind(element);
  }

public:
  int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const override final
  {
    return local_binary_integrand_->order(basis, *local_function_, param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalBasisType& basis,
                const DomainType& point_in_reference_element,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << this->logging_id << ".evaluate(basis.size()=" << basis.size(param)
                << ", point_in_reference_element=[" << point_in_reference_element << "], param=" << param
                << "):" << std::endl;
    // prepare storage
    const auto size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    // evaluate
    local_binary_integrand_->evaluate(
        basis, *local_function_, point_in_reference_element, local_binary_integrand_result_, param);
    // extract result
    LOG_(debug) << "  local_binary_integrand_result_ = [" << local_binary_integrand_result_ << "]" << std::endl;
    assert(local_binary_integrand_result_.rows() >= size && "This must not happen!");
    assert(local_binary_integrand_result_.cols() >= 1 && "This must not happen!");
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] = local_binary_integrand_result_[ii][0];
  } // ... evaluate(...)

private:
  const XT::Functions::GridFunction<E, a_r, a_rC, AF> inducing_function_as_ansatz_basis_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, a_r, a_rC, AF>::LocalFunctionType> local_function_;
  std::unique_ptr<LocalBinaryElementIntegrandType> local_binary_integrand_;
  mutable DynamicMatrix<F> local_binary_integrand_result_;
}; // class LocalBinaryToUnaryElementIntegrand


/**
 * \sa LocalBinaryToUnaryElementIntegrand
 * \sa LocalUnaryIntersectionIntegrandInterface
 * \sa LocalBinaryIntersectionIntegrandInterface
 * \sa XT::Functions::GridFunction
 */
template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryToUnaryIntersectionIntegrand : public LocalUnaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F>
{
  using ThisType = LocalBinaryToUnaryIntersectionIntegrand<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using BaseType = LocalUnaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalBasisType;

  using LocalBinaryIntersectionIntegrandType =
      LocalBinaryIntersectionIntegrandInterface<I, t_r, t_rC, TF, F, a_r, a_rC, AF>;

  LocalBinaryToUnaryIntersectionIntegrand(
      const LocalBinaryIntersectionIntegrandType& local_binary_integrand,
      XT::Functions::GridFunction<E, a_r, a_rC, AF> inducing_function_as_ansatz_basis)
    : inducing_function_as_ansatz_basis_(inducing_function_as_ansatz_basis)
    , local_function_(inducing_function_as_ansatz_basis_.local_function())
    , local_binary_integrand_(local_binary_integrand.copy_as_binary_intersection_integrand())
  {}

  LocalBinaryToUnaryIntersectionIntegrand(const ThisType& other)
    : BaseType()
    , inducing_function_as_ansatz_basis_(other.inducing_function_as_ansatz_basis_)
    , local_function_(inducing_function_as_ansatz_basis_.local_function())
    , local_binary_integrand_(other.local_binary_integrand_->copy_as_binary_intersection_integrand())
  {}

  std::unique_ptr<BaseType> copy_as_unary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  bool inside() const override final
  {
    return local_binary_integrand_->inside();
  }

protected:
  void post_bind(const IntersectionType& intersctn) override final
  {
    local_function_->bind(this->inside() ? intersctn.inside() : intersctn.outside());
    local_binary_integrand_->bind(intersctn);
  }

public:
  int order(const LocalBasisType& test_basis, const XT::Common::Parameter& param = {}) const override final
  {
    return local_binary_integrand_->order(test_basis, *local_function_, param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalBasisType& test_basis,
                const DomainType& point_in_reference_Intersection,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const auto size = test_basis.size(param);
    if (result.size() < size)
      result.resize(size);
    // evaluate
    local_binary_integrand_->evaluate(
        test_basis, *local_function_, point_in_reference_Intersection, local_binary_integrand_result_, param);
    // extract result
    assert(local_binary_integrand_result_.rows() >= size && "This must not happen!");
    assert(local_binary_integrand_result_.cols() >= 1 && "This must not happen!");
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] = local_binary_integrand_result_[ii][0];
  } // ... evaluate(...)

private:
  const XT::Functions::GridFunction<E, a_r, a_rC, AF> inducing_function_as_ansatz_basis_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, a_r, a_rC, AF>::LocalFunctionType> local_function_;
  std::unique_ptr<LocalBinaryIntersectionIntegrandType> local_binary_integrand_;
  mutable DynamicMatrix<F> local_binary_integrand_result_;
}; // class LocalBinaryToUnaryIntersectionIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
