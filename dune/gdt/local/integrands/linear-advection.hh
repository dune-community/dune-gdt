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

#include <dune/gdt/print.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * Given a direction v, computes
 *   - phi * (v \cdot \nabla psi), if advection_in_divergence_form
 * and
 *   (v \cdot \nabla phi) * psi
 * else, for all combinations of ansatz basis functions phi and test basis functions psi.
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

  explicit LocalLinearAdvectionIntegrand(XT::Functions::GridFunction<E, d, 1, F> direction,
                                         bool advection_in_divergence_form = true,
                                         const std::string& logging_prefix = "",
                                         const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType(direction.parameter_type(),
               logging_prefix.empty() ? "LocalLinearAdvectionIntegrand" : logging_prefix,
               logging_state)
    , direction_(direction.copy_as_grid_function())
    , advection_in_divergence_form_(advection_in_divergence_form)
    , local_direction_(direction_->local_function())
  {
    LOG_(info) << "LocalLinearAdvectionIntegrand(this=" << this << ", direction=" << &direction << ")" << std::endl;
  }

  LocalLinearAdvectionIntegrand(const ThisType& other)
    : BaseType(other)
    , direction_(other.direction_->copy_as_grid_function())
    , advection_in_divergence_form_(other.advection_in_divergence_form_)
    , local_direction_(direction_->local_function())
  {}

  LocalLinearAdvectionIntegrand(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_element_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const ElementType& ele) override
  {
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
    LOG_(debug) << "evaluate(test_basis.size()=" << test_basis.size(param)
                << ", ansatz_basis.size()=" << ansatz_basis.size(param)
                << ",\n    point_in_{reference_element|physical_space}={" << print(point_in_reference_element) << "|"
                << print(this->element().geometry().global(point_in_reference_element)) << "}"
                << ",\n    param=" << print(param) << ")" << std::endl;
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
    if (advection_in_divergence_form_) {
      test_basis.jacobians(point_in_reference_element, test_basis_grads_, param);
      ansatz_basis.evaluate(point_in_reference_element, ansatz_basis_values_, param);
      LOG_(debug) << "  test_basis_grads_ = " << print(test_basis_grads_, {{"oneline", "true"}})
                  << "\n  ansatz_basis_values_ = " << print(ansatz_basis_values_, {{"oneline", "true"}})
                  << "\n  direction = " << direction << std::endl;
      // compute integrand
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += -1.0 * ansatz_basis_values_[jj] * (direction * test_basis_grads_[ii][0]);
    } else {
      test_basis.evaluate(point_in_reference_element, test_basis_values_, param);
      ansatz_basis.jacobians(point_in_reference_element, ansatz_basis_grads_, param);
      LOG_(debug) << "  test_basis_values_ = " << print(test_basis_values_, {{"oneline", "true"}})
                  << "\n  ansatz_basis_grads_ = " << print(ansatz_basis_grads_, {{"oneline", "true"}})
                  << "\n  direction = " << direction << std::endl;
      // compute integrand
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += (direction * ansatz_basis_grads_[jj][0]) * test_basis_values_[ii];
    }
    LOG_(debug) << "  result = " << print(result, {{"oneline", "true"}}) << std::endl;
  } // ... evaluate(...)

private:
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, d, 1, F>> direction_;
  const bool advection_in_divergence_form_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, 1, F>::LocalFunctionType> local_direction_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
}; // class LocalLinearAdvectionIntegrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_HH
