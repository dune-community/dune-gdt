// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ABS_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ABS_HH

#include <cmath>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class E, size_t r = 1, size_t rC = 1, class R = double, class F = double>
class LocalElementAbsIntegrand : public LocalUnaryElementIntegrandInterface<E, r, rC, R, F>
{
  static_assert(rC == 1, "");
  using ThisType = LocalElementAbsIntegrand;
  using BaseType = LocalUnaryElementIntegrandInterface<E, r, rC, R, F>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::LocalBasisType;

  std::unique_ptr<BaseType> copy_as_unary_element_integrand() const
  {
    return std::make_unique<ThisType>();
  }

  int order(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const
  {
    return basis.order(param);
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
    basis.evaluate(point_in_reference_element, basis_values_, param);
    // compute integrand
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] = basis_values_[ii].two_norm();
  } // ... evaluate(...)

private:
  mutable std::vector<typename LocalBasisType::RangeType> basis_values_;
}; // class LocalElementAbsIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ABS_HH
