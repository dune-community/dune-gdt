// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_UPWIND_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_UPWIND_HH

#include <dune/xt/functions/grid-function.hh>

#include <dune/gdt/print.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace LocalLinearAdvectionUpwindIntegrands {


template <class E, class F = double>
class Volume : public LocalBinaryElementIntegrandInterface<E, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalBinaryElementIntegrandInterface<E, 1, 1, F, F, 1, 1, F>;
  using ThisType = Volume;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  explicit Volume(XT::Functions::GridFunction<E, d, 1, F> direction = FieldVector<F, d>(1.),
                  const std::string& logging_prefix = "")
    : BaseType(direction.parameter_type(),
               logging_prefix.empty() ? "gdt" : "gdt.locallinearadvectionupwindvolumeintegrand",
               logging_prefix.empty() ? "LocalLinearAdvectionUpwindIntegrands::Volume" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , direction_(direction)
    , local_direction_(direction_.local_function())
  {
    LOG_(info) << this->logging_id << "(direction=" << &direction << ")" << std::endl;
  }

  Volume(const ThisType& other)
    : BaseType(other.parameter_type())
    , direction_(other.direction_)
    , local_direction_(direction_.local_function())
  {}

  Volume(ThisType&& source) = default;

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
    LOG_(debug) << "evaluate(test_basis.size()=" << test_basis.size(param)
                << ", ansatz_basis.size()=" << ansatz_basis.size(param)
                << ", point_in_reference_element=" << print(point_in_reference_element) << ", param=" << param << ")"
                << std::endl;
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // evaluate
    test_basis.evaluate(point_in_reference_element, test_basis_values_, param);
    ansatz_basis.jacobians(point_in_reference_element, ansatz_basis_grads_, param);
    const auto direction = local_direction_->evaluate(point_in_reference_element, param);
    LOG_(debug) << "  test_basis_values_ = " << print(test_basis_values_)
                << "\n  ansatz_basis_grads_ = " << print(ansatz_basis_grads_) << "\n  direction = " << direction
                << std::endl;
    // compute integrand
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += (direction * ansatz_basis_grads_[ii][0]) * test_basis_values_[jj];
    LOG_(debug) << "  result = " << result << std::endl;
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, d, 1, F> direction_;
  std::unique_ptr<typename XT::Functions::GridFunction<E, d, 1, F>::LocalFunctionType> local_direction_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class Volume


/**
 * \note Makes sense on outflow intersections, i.e. where direction*normal < 0
 * \sa LocalCouplingIntersectionRestrictedIntegralBilinearForm
 */
template <class I>
class InnerCoupling : public LocalQuaternaryIntersectionIntegrandInterface<I>
{
  using ThisType = InnerCoupling;
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I>;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  InnerCoupling(XT::Functions::GridFunction<E, d> direction, const std::string& logging_prefix = "")
    : BaseType(direction.parameter_type(),
               logging_prefix.empty() ? "gdt" : "gdt.locallinearadvectionupwindinnercouplingintegrand",
               logging_prefix.empty() ? "LocalLinearAdvectionUpwindIntegrands::InnerCoupling" : logging_prefix,
               /*logging_disabled=*/logging_prefix.empty())
    , direction_(direction)
    , local_direction_in_(direction_.local_function())
  {
    LOG_(info) << this->logging_id << "(direction=" << &direction << ")" << std::endl;
  }

  InnerCoupling(const ThisType& other)
    : BaseType(other.parameter_type())
    , direction_(other.direction_)
    , local_direction_in_(direction_.local_function())
  {}

  InnerCoupling(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) override final
  {
    DUNE_THROW_IF(
        !intrsctn.neighbor(), Exceptions::integrand_error, "This integrand cannot be used on a boundary intersection!");
    const auto inside_element = intrsctn.inside();
    local_direction_in_->bind(inside_element);
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& /*test_basis_outside*/,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_direction_in_->order(param) + test_basis_inside.order(param)
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
    LOG_(debug) << "evaluate(test_basis_{inside,outside}.size()={" << test_basis_inside.size(param) << ","
                << test_basis_outside.size(param) << "}, ansatz_basis_{inside,outside}.size()={"
                << ansatz_basis_inside.size(param) << "," << ansatz_basis_outside.size(param)
                << "}, point_in_reference_intersection=" << print(point_in_reference_intersection)
                << ", param=" << param << ")" << std::endl;
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
    // ... basis functions  ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_in_values_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_in_values_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    // ... data functions, ...
    const auto direction = local_direction_in_->evaluate(point_in_inside_reference_element, param);
    LOG_(debug) << "  test_basis_in_values_ = " << print(test_basis_in_values_)
                << "\n  ansatz_basis_in_values_ = " << print(ansatz_basis_in_values_)
                << "\n  ansatz_basis_out_values_ = " << print(ansatz_basis_out_values_)
                << "\n  normal = " << print(normal) << "\n  direction = " << direction << std::endl;
    // ... and finally compute the integrand.
    const size_t rows_in = test_basis_inside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_in_in[ii][jj] += -1.0 * (direction * normal) * ansatz_basis_in_values_[jj] * test_basis_in_values_[ii];
      for (size_t jj = 0; jj < cols_out; ++jj)
        result_out_in[ii][jj] += (direction * normal) * ansatz_basis_out_values_[jj] * test_basis_in_values_[ii];
    }
    // nothing to do for test_basis_outside ...
    LOG_(debug) << "  result_in_in = " << result_in_in << std::endl;
    LOG_(debug) << "  result_out_in = " << result_out_in << std::endl;
  } // ... evaluate(...)

private:
  XT::Functions::GridFunction<E, d> direction_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d>::LocalFunctionType> local_direction_in_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; // InnerCoupling


/**
 * \note Makes sense on outflow intersections, i.e. where direction*normal < 0
 * \sa LocalIntersectionRestrictedIntegralBilinearForm
 */
template <class I>
class DirichletCoupling
  : public LocalUnaryIntersectionIntegrandInterface<I>
  , public LocalBinaryIntersectionIntegrandInterface<I>
{
  using ThisType = DirichletCoupling;
  using BaseUnaryType = LocalUnaryIntersectionIntegrandInterface<I>;
  using BaseBinaryType = LocalBinaryIntersectionIntegrandInterface<I>;

public:
  using BaseBinaryType::d;
  using typename BaseBinaryType::DomainType;
  using typename BaseBinaryType::E;
  using typename BaseBinaryType::F;
  using typename BaseBinaryType::IntersectionType;
  using typename BaseBinaryType::LocalAnsatzBasisType;
  using typename BaseBinaryType::LocalTestBasisType;

  /**
   * \note dirichlet_data is only required if used as a unary integrand, i.e. for the right hand side
   */
  DirichletCoupling(XT::Functions::GridFunction<E, d> direction,
                    XT::Functions::GridFunction<E> dirichlet_data = 0.,
                    const std::string& logging_prefix = "")
    : BaseUnaryType(direction.parameter_type() + dirichlet_data.parameter_type(),
                    logging_prefix.empty() ? "gdt" : "gdt.locallinearadvectionupwinddirichletcouplingintegrand",
                    logging_prefix.empty() ? "LocalLinearAdvectionUpwindIntegrands::DirichletCoupling" : logging_prefix,
                    /*logging_disabled=*/logging_prefix.empty())
    , BaseBinaryType(direction.parameter_type() + dirichlet_data.parameter_type(),
                     logging_prefix.empty() ? "gdt" : "gdt.locallinearadvectionupwinddirichletcouplingintegrand",
                     logging_prefix.empty() ? "LocalLinearAdvectionUpwindIntegrands::DirichletCoupling"
                                            : logging_prefix,
                     /*logging_disabled=*/logging_prefix.empty())
    , direction_(direction)
    , dirichlet_data_(dirichlet_data)
    , local_direction_(direction_.local_function())
    , local_dirichlet_data_(dirichlet_data_.local_function())
  {
    LOG_(info) << this->logging_id << "(direction=" << &direction << ", dirichlet_data=" << &dirichlet_data << ")"
               << std::endl;
  }

  DirichletCoupling(const ThisType& other)
    : BaseUnaryType(other)
    , BaseBinaryType(other)
    , direction_(other.direction_)
    , dirichlet_data_(other.dirichlet_data_)
    , local_direction_(direction_.local_function())
    , local_dirichlet_data_(dirichlet_data_.local_function())
  {}

  DirichletCoupling(ThisType&& source) = default;

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    const auto inside_element = intersection.inside();
    local_direction_->bind(inside_element);
    local_dirichlet_data_->bind(inside_element);
  }

public:
  bool inside() const override final
  {
    return true; // We expect the bases to be bound to the inside (see evaluate and post_bind).
  }

  /// \name Required by LocalUnaryIntersectionIntegrandInterface.
  /// \{

  std::unique_ptr<BaseUnaryType> copy_as_unary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalTestBasisType& test_basis, const XT::Common::Parameter& param = {}) const override final
  {
    return local_dirichlet_data_->order(param) + local_direction_->order(param) + test_basis.order(param);
  }

  void evaluate(const LocalTestBasisType& test_basis,
                const DomainType& point_in_reference_intersection,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "evaluate(test_basis.size()=" << test_basis.size(param)
                << ", point_in_reference_intersection=" << print(point_in_reference_intersection) << ", param=" << param
                << ")" << std::endl;
    // Prepare sotrage, ...
    BaseUnaryType::ensure_size_and_clear_results(test_basis, result, param);
    // evaluate ...
    const auto point_in_inside_reference_element =
        BaseUnaryType::intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = BaseUnaryType::intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    // ... data functions, ...
    const auto direction = local_direction_->evaluate(point_in_inside_reference_element, param);
    const auto dirichlet_data = local_dirichlet_data_->evaluate(point_in_inside_reference_element, param);
    LOG_(debug) << "  test_basis_values_ = " << print(test_basis_values_) << "\n  normal = " << print(normal)
                << "\n  direction = " << direction << "\n  dirichlet_data = " << dirichlet_data << std::endl;
    // ... and finally compute the integrand.
    const size_t size = test_basis.size(param);
    for (size_t jj = 0; jj < size; ++jj)
      result[jj] += -1.0 * (direction * normal) * dirichlet_data * test_basis_values_[jj];
    LOG_(debug) << "  result = " << result << std::endl;
  } // ... evaluate(...)

  /// \}
  /// \name Required by LocalBinaryIntersectionIntegrandInterface.
  /// \{

  std::unique_ptr<BaseBinaryType> copy_as_binary_intersection_integrand() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_direction_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "evaluate(test_basis.size()=" << test_basis.size(param)
                << ", ansatz_basis.size()=" << ansatz_basis.size(param)
                << ", point_in_reference_intersection=" << print(point_in_reference_intersection) << ", param=" << param
                << ")" << std::endl;
    // Prepare sotrage, ...
    BaseBinaryType::ensure_size_and_clear_results(test_basis, ansatz_basis, result, param);
    // evaluate ...
    const auto point_in_inside_reference_element =
        BaseBinaryType::intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = BaseBinaryType::intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    ansatz_basis.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    // ... data functions, ...
    const auto direction = local_direction_->evaluate(point_in_inside_reference_element, param);
    LOG_(debug) << "  test_basis_values_ = " << print(test_basis_values_)
                << "\n  ansatz_basis_values_ = " << print(ansatz_basis_values_) << "\n  normal = " << print(normal)
                << "\n  direction = " << direction << std::endl;
    // ... and finally compute the integrand.
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        result[ii][jj] += -1.0 * (direction * normal) * ansatz_basis_values_[jj] * test_basis_values_[ii];
    LOG_(debug) << "  result = " << result << std::endl;
  } // ... evaluate(...)

  /// \}

private:
  XT::Functions::GridFunction<E, d> direction_;
  XT::Functions::GridFunction<E> dirichlet_data_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d>::LocalFunctionType> local_direction_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E>::LocalFunctionType> local_dirichlet_data_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
}; // class DirichletCoupling


} // namespace LocalLinearAdvectionUpwindIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LINEAR_ADVECTION_UPWIND_HH