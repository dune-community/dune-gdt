// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH

#include <dune/xt/functions/interfaces/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given a function f and a binary element integrand bi(\cdot, cdot), models the unary element integrand bi(f, \cdot).
 *
 * See also LocalUnaryElementIntegrandInterface for a description of the template arguments.
 *
 * \sa local_binary_to_unary_element_integrand
 * \sa LocalUnaryElementIntegrandInterface
 * \sa LocalBinaryElementIntegrandInterface
 * \sa XT::Functions::GridFunctionInterface
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TF = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AF = TF>
class LocalBinaryToUnaryElementIntegrand : public LocalUnaryElementIntegrandInterface<E, a_r, a_rC, AF, F>
{
  using ThisType = LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;
  using BaseType = LocalUnaryElementIntegrandInterface<E, a_r, a_rC, AF, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalBasisType;

  using LocalizableFunctionType = XT::Functions::GridFunctionInterface<E, t_r, t_rC, TF>;
  using LocalBinaryElementIntegrandType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>;

  LocalBinaryToUnaryElementIntegrand(const LocalizableFunctionType& inducing_function_as_test_basis,
                                     const LocalBinaryElementIntegrandType& local_binary_integrand)
    : inducing_function_as_test_basis_(inducing_function_as_test_basis)
    , local_function_(inducing_function_as_test_basis_.local_function())
    , local_binary_integrand_(local_binary_integrand.copy())
  {
  }

  LocalBinaryToUnaryElementIntegrand(const ThisType& other)
    : BaseType()
    , inducing_function_as_test_basis_(other.inducing_function_as_test_basis_)
    , local_function_(inducing_function_as_test_basis_.local_function())
    , local_binary_integrand_(other.local_binary_integrand_->copy())
  {
  }

  std::unique_ptr<BaseType> copy() const override final
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
    return local_binary_integrand_->order(*local_function_, basis, param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalBasisType& basis,
                const DomainType& point_in_reference_element,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const auto size = basis.size(param);
    if (result.size() < size)
      result.resize(size);
    // evaluate
    local_binary_integrand_->evaluate(
        *local_function_, basis, point_in_reference_element, local_binary_integrand_result_, param);
    // extract result
    result = local_binary_integrand_result_[0];
  } // ... evaluate(...)

private:
  const LocalizableFunctionType& inducing_function_as_test_basis_;
  std::unique_ptr<typename LocalizableFunctionType::LocalFunctionType> local_function_;
  std::unique_ptr<LocalBinaryElementIntegrandType> local_binary_integrand_;
  mutable DynamicMatrix<F> local_binary_integrand_result_;
}; // class LocalBinaryToUnaryElementIntegrand


/**
 * \sa LocalBinaryToUnaryElementIntegrand
 */
template <class E, size_t t_r, size_t t_rC, class TF, class F, size_t a_r, size_t a_rC, class AF>
LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TF, F, a_r, a_rC, AF> local_binary_to_unary_element_integrand(
    const XT::Functions::GridFunctionInterface<E, t_r, t_rC, TF>& inducing_function_as_test_basis,
    const LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TF, F, a_r, a_rC, AF>& local_binary_element_integrand)
{
  return LocalBinaryToUnaryElementIntegrand<E, t_r, t_rC, TF, F, a_r, a_rC, AF>(inducing_function_as_test_basis,
                                                                                local_binary_element_integrand);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_CONVERSION_HH
