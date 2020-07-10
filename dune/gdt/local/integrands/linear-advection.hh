// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/functions/grid-function.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given a direction v, computes `- (v phi) * \nabla psi` for all combinations of ansatz basis functions phi and test
 * basis functions psi.
 */
template <class E, class F = double>
class LocalLinearAdvectionIntegrand : public LocalBinaryElementIntegrandInterface<E, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, 1, 1, F, F, 1, 1, F>;
  using ThisType = LocalLinearAdvectionIntegrand;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  explicit LocalLinearAdvectionIntegrand(XT::Functions::GridFunction<E, d, 1, F> direction = FieldVector<F, d>(1.),
                                         const std::string& logging_prefix = "")
    : BaseType(direction.parameter_type(),
               logging_prefix.empty() ? "gdt" : "gdt.locallinearadvectionintegrand",
               logging_prefix.empty() ? "LocalLinearAdvectionIntegrand" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , direction_(direction)
    , local_direction_(direction_.local_function())
  {}

  LocalLinearAdvectionIntegrand(const ThisType& other)
    : BaseType(other)
    , direction_(other.direction_)
    , local_direction_(direction_.local_function())
  {}

  LocalLinearAdvectionIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override
  {
#ifndef NDEBUG
    if (!ele.geometry().affine())
      std::cerr << "Warning: integration order has to be increased for non-affine geometries!" << std::endl;
#endif
    local_direction_->bind(ele);
  }

public:
  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_direction_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

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
    result *= 0;
    // evaluate
    test_basis.jacobians(point_in_reference_element, test_basis_grads_, param);
    ansatz_basis.evaluate(point_in_reference_element, ansatz_basis_values_, param);
    const auto direction = local_direction_->evaluate(point_in_reference_element, param);
    // compute integrand
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += -1.0 * ansatz_basis_values_[jj] * (direction * test_basis_grads_[ii][0]);
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, d, 1, F> direction_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, 1, F>::LocalFunctionType> local_direction_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // class LocalLinearAdvectionIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_HH
