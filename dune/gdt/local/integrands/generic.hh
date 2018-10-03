// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

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
  using typename BaseType::LocalBasisType;
  using typename BaseType::DomainType;

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
  {
  }

  GenericLocalUnaryElementIntegrand(const ThisType& other)
    : BaseType(other.parameter_type())
    , order_(other.order_)
    , evaluate_(other.evaluate_)
  {
  }

  std::unique_ptr<BaseType> copy() const
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


} // namespace GDT
} // namespace Dun

#endif // DUNE_GDT_LOCAL_INTEGRANDS_GENERIC_HH
