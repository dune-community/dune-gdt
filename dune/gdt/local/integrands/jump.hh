// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_JUMP_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_JUMP_HH

#include <functional>

#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

#include "interfaces.hh"
#include "ipdg.hh"

namespace Dune {
namespace GDT {
namespace LocalJumpIntegrands {


template <class I, size_t r = 1>
class Inner : public LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, double, double, r, 1>
{
public:
  using ThisType = Inner;
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, r, 1, double, double, r, 1>;

  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  Inner(const std::function<double(const I&)>& intersection_diameter =
            LocalIPDGIntegrands::internal::default_intersection_diameter<I>(),
        const std::string& logging_prefix = "",
        const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({},
               logging_prefix.empty() ? "LocalJumpIntegrands::Inner" : logging_prefix,
               logging_state)
    , intersection_diameter_(intersection_diameter)
  {}

  Inner(const ThisType& other)
    : BaseType(other)
    , intersection_diameter_(other.intersection_diameter_)
  {}

  Inner(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    DUNE_THROW_IF(
        !intrsctn.neighbor(), Exceptions::integrand_error, "This integrand cannot be used on a boundary intersection!");
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(test_basis_inside.order(param), test_basis_outside.order(param))
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
    // Prepare sotrage, ...
    this->ensure_size_and_clear_results(test_basis_inside,
                                        ansatz_basis_inside,
                                        test_basis_outside,
                                        ansatz_basis_outside,
                                        result_in_in,
                                        result_in_out,
                                        result_out_in,
                                        result_out_out,
                                        param);
    // evaluate ...
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    const auto point_in_outside_reference_element =
        this->intersection().geometryInOutside().global(point_in_reference_intersection);
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_in_values_, param);
    test_basis_outside.evaluate(point_in_outside_reference_element, test_basis_out_values_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_in_values_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    // ... and diameter
    const auto h = intersection_diameter_(this->intersection());
    // and finally compute the integrand.
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    if constexpr (r == 1) {
      for (size_t ii = 0; ii < rows_in; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_in_in[ii][jj] += h * ansatz_basis_in_values_[jj] * test_basis_in_values_[ii];
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_in_out[ii][jj] += -1.0 * h * ansatz_basis_out_values_[jj] * test_basis_in_values_[ii];
      }
      for (size_t ii = 0; ii < rows_out; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_out_in[ii][jj] += -1.0 * h * ansatz_basis_in_values_[jj] * test_basis_out_values_[ii];
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_out_out[ii][jj] += h * ansatz_basis_out_values_[jj] * test_basis_out_values_[ii];
      }
    } else {
      for (size_t ii = 0; ii < rows_in; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_in_in[ii][jj] += h * (ansatz_basis_in_values_[jj] * normal) * (test_basis_in_values_[ii] * normal);
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_in_out[ii][jj] +=
              -1.0 * h * (ansatz_basis_out_values_[jj] * normal) * (test_basis_in_values_[ii] * normal);
      }
      for (size_t ii = 0; ii < rows_out; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_out_in[ii][jj] +=
              -1.0 * h * (ansatz_basis_in_values_[jj] * normal) * (test_basis_out_values_[ii] * normal);
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_out_out[ii][jj] += h * (ansatz_basis_out_values_[jj] * normal) * (test_basis_out_values_[ii] * normal);
      }
    }
  } // ... evaluate(...)

private:
  const std::function<double(const I&)> intersection_diameter_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; // Inner


template <class I, size_t r = 1>
class Boundary : public LocalBinaryIntersectionIntegrandInterface<I, r, 1, double, double, r, 1>
{
public:
  using ThisType = Boundary;
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I, r, 1, double, double, r, 1>;

  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  Boundary(const std::function<double(const I&)>& intersection_diameter =
               LocalIPDGIntegrands::internal::default_intersection_diameter<I>(),
           const std::string& logging_prefix = "",
           const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({},
               logging_prefix.empty() ? "LocalJumpIntegrands::Boundary" : logging_prefix,
               logging_state)
    , intersection_diameter_(intersection_diameter)
  {}

  Boundary(const ThisType& other)
    : BaseType(other)
    , intersection_diameter_(other.intersection_diameter_)
  {}

  Boundary(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  bool inside() const override final
  {
    return true; // We expect the bases to be bound to the inside (see evaluate).
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return test_basis.order(param) + ansatz_basis.order(param);
  }

  using BaseType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    // Prepare sotrage, ...
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // evaluate ...
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions ...
    test_basis.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    ansatz_basis.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    //  ... and diameter ...
    const auto h = intersection_diameter_(this->intersection());
    // and finally compute integrand.
    if constexpr (r == 1)
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += h * ansatz_basis_values_[jj] * test_basis_values_[ii];
    else
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += h * (ansatz_basis_values_[jj] * normal) * (test_basis_values_[ii] * normal);
  } // ... evaluate(...)

private:
  const std::function<double(const I&)> intersection_diameter_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // Boundary


} // namespace LocalJumpIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_JUMP_HH
