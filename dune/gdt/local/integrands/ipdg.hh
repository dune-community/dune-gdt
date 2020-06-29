// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_IPDG_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_IPDG_HH

#include <functional>

#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/intersection.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace LocalIPDGIntegrands {
namespace internal {


template <class Intersection>
static std::function<double(const Intersection&)> default_inner_intersection_diameter()
{
  return [](const Intersection& intersection) {
    if (Intersection::dimension == 1) {
      if (intersection.neighbor())
        return 0.5 * (XT::Grid::diameter(intersection.inside()) + XT::Grid::diameter(intersection.outside()));
      else
        return XT::Grid::diameter(intersection.inside());
    } else
      return XT::Grid::diameter(intersection);
  };
} // ... default_inner_intersection_diameter(...)


template <class Intersection>
static std::function<double(const Intersection&)> default_boundary_intersection_diameter()
{
  return [](const Intersection& intersection) {
    if (Intersection::dimension == 1)
      return XT::Grid::diameter(intersection.inside());
    else
      return XT::Grid::diameter(intersection);
  };
} // ... default_boundary_intersection_diameter(...)


} // namespace internal


template <class I>
class InnerPenalty : public LocalQuaternaryIntersectionIntegrandInterface<I>
{
  using ThisType = InnerPenalty;
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  InnerPenalty(
      const double& penalty,
      XT::Functions::GridFunction<E, d, d> weight_function = 1.,
      const std::function<double(const I&)>& intersection_diameter = internal::default_inner_intersection_diameter<I>())
    : BaseType(weight_function.parameter_type())
    , penalty_(penalty)
    , weight_(weight_function)
    , intersection_diameter_(intersection_diameter)
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  InnerPenalty(const ThisType& other)
    : BaseType(other.parameter_type())
    , penalty_(other.penalty_)
    , weight_(other.weight_)
    , intersection_diameter_(other.intersection_diameter_)
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  InnerPenalty(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    DUNE_THROW_IF(
        !intrsctn.neighbor(), Exceptions::integrand_error, "This integrand cannot be used on a boundary intersection!");
    local_weight_in_->bind(intrsctn.inside());
    local_weight_out_->bind(intrsctn.outside());
  }

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
    // ... and data functions, ...
    const auto weight_in = local_weight_in_->evaluate(point_in_inside_reference_element, param);
    const auto weight_out = local_weight_out_->evaluate(point_in_outside_reference_element, param);
    // compute the weighted penalty ...
    const double delta_plus = normal * (weight_out * normal);
    const double delta_minus = normal * (weight_in * normal);
    const auto weight = (delta_plus * delta_minus) / (delta_plus + delta_minus); // half harmonic average
    const auto h = intersection_diameter_(this->intersection());
    const auto penalty = (penalty_ * weight) / h;
    // and finally compute the integrand.
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_in_in[ii][jj] += penalty * ansatz_basis_in_values_[jj] * test_basis_in_values_[ii];
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_in_out[ii][jj] += -1.0 * penalty * ansatz_basis_out_values_[jj] * test_basis_in_values_[ii];
    }
    for (size_t ii = 0; ii < rows_out; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_out_in[ii][jj] += -1.0 * penalty * ansatz_basis_in_values_[jj] * test_basis_out_values_[ii];
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_out_out[ii][jj] += penalty * ansatz_basis_out_values_[jj] * test_basis_out_values_[ii];
    }
  } // ... evaluate(...)

private:
  const double penalty_;
  XT::Functions::GridFunction<E, d, d> weight_;
  const std::function<double(const I&)> intersection_diameter_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, d>::LocalFunctionType> local_weight_in_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, d>::LocalFunctionType> local_weight_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; // InnerPenalty


template <class I>
class BoundaryPenalty : public LocalBinaryIntersectionIntegrandInterface<I>
{
  using ThisType = BoundaryPenalty;
  using BaseType = LocalBinaryIntersectionIntegrandInterface<I>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  BoundaryPenalty(const double& penalty,
                  XT::Functions::GridFunction<E, d, d> weight_function = 1.,
                  const std::function<double(const I&)>& intersection_diameter =
                      internal::default_boundary_intersection_diameter<I>())
    : BaseType()
    , penalty_(penalty)
    , weight_(weight_function)
    , intersection_diameter_(intersection_diameter)
    , local_weight_(weight_.local_function())
  {}

  BoundaryPenalty(const ThisType& other)
    : BaseType(other)
    , penalty_(other.penalty_)
    , weight_(other.weight_)
    , intersection_diameter_(other.intersection_diameter_)
    , local_weight_(weight_.local_function())
  {}

  BoundaryPenalty(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    local_weight_->bind(intersection.inside());
  }

public:
  bool inside() const override final
  {
    return true; // We expect the bases to be bound to the inside (see evaluate).
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
    // ... and data functions, ....
    const auto weight = local_weight_->evaluate(point_in_inside_reference_element, param);
    // compute the weighted penalty ...
    const auto h = intersection_diameter_(this->intersection());
    const auto penalty = (penalty_ * (normal * (weight * normal))) / h;
    // and finally compute integrand.
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += penalty * ansatz_basis_values_[jj] * test_basis_values_[ii];
  } // ... evaluate(...)

private:
  const double penalty_;
  XT::Functions::GridFunction<E, d, d> weight_;
  const std::function<double(const I&)> intersection_diameter_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, d>::LocalFunctionType> local_weight_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // BoundaryPenalty


} // namespace LocalIPDGIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_IPDG_HH
