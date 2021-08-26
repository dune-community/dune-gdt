// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_IPDG_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_IPDG_HH

#include <dune/xt/functions/grid-function.hh>

#include "interfaces.hh"
#include "ipdg.hh"

namespace Dune {
namespace GDT {
namespace LocalLaplaceIPDGIntegrands {


/**
 * \note The role of symmetry_prefactor:
 *       * -1 => NIPDG
 *       *  0 => IIPDG
 *       *  1 => SIPDG
 * \note The role of the weight:
 *       * symmetry_prefactor = 1 && weight_function = 1 => SIPDG
 *       * symmetry_prefactor = 1 && weight_function = diffusion => SWIPDG
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

  InnerCoupling(const double& symmetry_prefactor,
                XT::Functions::GridFunction<E, d, d> diffusion,
                XT::Functions::GridFunction<E, d, d> weight_function = {1.},
                const std::string& logging_prefix = "",
                const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType(diffusion.parameter_type() + weight_function.parameter_type(),
               logging_prefix.empty() ? "LocalLaplaceIPDGIntegrands::InnerCoupling" : logging_prefix,
               logging_state)
    , symmetry_prefactor_(symmetry_prefactor)
    , diffusion_(diffusion.copy_as_grid_function())
    , weight_(weight_function.copy_as_grid_function())
    , local_diffusion_in_(diffusion_->local_function())
    , local_diffusion_out_(diffusion_->local_function())
    , local_weight_in_(weight_->local_function())
    , local_weight_out_(weight_->local_function())
  {}

  InnerCoupling(const ThisType& other)
    : BaseType(other)
    , symmetry_prefactor_(other.symmetry_prefactor_)
    , diffusion_(other.diffusion_->copy_as_grid_function())
    , weight_(other.weight_->copy_as_grid_function())
    , local_diffusion_in_(diffusion_->local_function())
    , local_diffusion_out_(diffusion_->local_function())
    , local_weight_in_(weight_->local_function())
    , local_weight_out_(weight_->local_function())
  {}

  InnerCoupling(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy_as_quaternary_intersection_integrand() const final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intrsctn) final
  {
    DUNE_THROW_IF(
        !intrsctn.neighbor(), Exceptions::integrand_error, "This integrand cannot be used on a boundary intersection!");
    const auto inside_element = intrsctn.inside();
    const auto outside_element = intrsctn.outside();
    local_diffusion_in_->bind(inside_element);
    local_weight_in_->bind(inside_element);
    local_diffusion_out_->bind(outside_element);
    local_weight_out_->bind(outside_element);
  } // ... post_bind(...)

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const final
  {
    return std::max(local_diffusion_in_->order(param), local_diffusion_out_->order(param))
           + std::max(local_weight_in_->order(), local_weight_out_->order(param))
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
                const XT::Common::Parameter& param = {}) const final
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
    // ... basis functions and ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_in_values_, param);
    test_basis_inside.jacobians(point_in_inside_reference_element, test_basis_in_grads_, param);
    test_basis_outside.evaluate(point_in_outside_reference_element, test_basis_out_values_, param);
    test_basis_outside.jacobians(point_in_outside_reference_element, test_basis_out_grads_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_in_values_, param);
    ansatz_basis_inside.jacobians(point_in_inside_reference_element, ansatz_basis_in_grads_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    ansatz_basis_outside.jacobians(point_in_outside_reference_element, ansatz_basis_out_grads_, param);
    // ... data functions, ...
    const auto diffusion_in = local_diffusion_in_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_out = local_diffusion_out_->evaluate(point_in_outside_reference_element, param);
    const auto weight_in = local_weight_in_->evaluate(point_in_inside_reference_element, param);
    const auto weight_out = local_weight_out_->evaluate(point_in_outside_reference_element, param);
    // compute the weighted mean ...
    const auto delta_plus = normal * (weight_out * normal);
    const auto delta_minus = normal * (weight_in * normal);
    const auto weight_minus = delta_plus / (delta_plus + delta_minus);
    const auto weight_plus = delta_minus / (delta_plus + delta_minus);
    // ... and finally compute the integrand.
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj) {
        result_in_in[ii][jj] +=
            -1.0 * weight_minus * ((diffusion_in * ansatz_basis_in_grads_[jj][0]) * normal) * test_basis_in_values_[ii];
        result_in_in[ii][jj] += -1.0 * symmetry_prefactor_ * weight_minus * ansatz_basis_in_values_[jj]
                                * ((diffusion_in * test_basis_in_grads_[ii][0]) * normal);
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        result_in_out[ii][jj] += -1.0 * weight_plus * ((diffusion_out * ansatz_basis_out_grads_[jj][0]) * normal)
                                 * test_basis_in_values_[ii];
        result_in_out[ii][jj] += symmetry_prefactor_ * weight_minus * ansatz_basis_out_values_[jj]
                                 * ((diffusion_in * test_basis_in_grads_[ii][0]) * normal);
      }
    }
    for (size_t ii = 0; ii < rows_out; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj) {
        result_out_in[ii][jj] +=
            weight_minus * ((diffusion_in * ansatz_basis_in_grads_[jj][0]) * normal) * test_basis_out_values_[ii];
        result_out_in[ii][jj] += -1.0 * symmetry_prefactor_ * weight_plus * ansatz_basis_in_values_[jj]
                                 * ((diffusion_out * test_basis_out_grads_[ii][0]) * normal);
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        result_out_out[ii][jj] +=
            weight_plus * ((diffusion_out * ansatz_basis_out_grads_[jj][0]) * normal) * test_basis_out_values_[ii];
        result_out_out[ii][jj] += symmetry_prefactor_ * weight_plus * ansatz_basis_out_values_[jj]
                                  * ((diffusion_out * test_basis_out_grads_[ii][0]) * normal);
      }
    }
  } // ... evaluate(...)

private:
  const double symmetry_prefactor_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, d, d>> diffusion_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, d, d>> weight_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_diffusion_in_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_diffusion_out_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_weight_in_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_weight_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_in_grads_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_out_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_in_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_out_grads_;
}; // InnerCoupling


/**
 * \note The role of symmetry_prefactor:
 *       * -1 => NIPDG
 *       *  0 => IIPDG
 *       *  1 => SIPDG
 * \note The role of the weight:
 *       * symmetry_prefactor = 1 && weight_function = 1 => SIPDG
 *       * symmetry_prefactor = 1 && weight_function = diffusion => SWIPDG
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
  DirichletCoupling(const double& symmetry_prefactor,
                    XT::Functions::GridFunction<E, d, d> diffusion,
                    XT::Functions::GridFunction<E> dirichlet_data = 0.,
                    const std::string& logging_prefix = "",
                    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseUnaryType(diffusion.parameter_type() + dirichlet_data.parameter_type(),
                    logging_prefix.empty() ? "LocalLaplaceIPDGIntegrands::DirichletCoupling" : logging_prefix,
                    logging_state)
    , BaseBinaryType(diffusion.parameter_type() + dirichlet_data.parameter_type(),
                     logging_prefix.empty() ? "LocalLaplaceIPDGIntegrands::DirichletCoupling" : logging_prefix,
                     logging_state)
    , symmetry_prefactor_(symmetry_prefactor)
    , diffusion_(diffusion.copy_as_grid_function())
    , dirichlet_data_(dirichlet_data.copy_as_grid_function())
    , local_diffusion_(diffusion_->local_function())
    , local_dirichlet_data_(dirichlet_data_->local_function())
  {}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
  DirichletCoupling(const ThisType& other)
    : BaseUnaryType(other)
    , BaseBinaryType(other)
    , symmetry_prefactor_(other.symmetry_prefactor_)
    , diffusion_(other.diffusion_->copy_as_grid_function())
    , dirichlet_data_(other.dirichlet_data_->copy_as_grid_function())
    , local_diffusion_(diffusion_->local_function())
    , local_dirichlet_data_(dirichlet_data_->local_function())
  {}
#pragma GCC diagnostic pop


  DirichletCoupling(ThisType&& source) = default;

protected:
  void post_bind(const IntersectionType& intersection) final
  {
    const auto inside_element = intersection.inside();
    local_diffusion_->bind(inside_element);
    local_dirichlet_data_->bind(inside_element);
  }

public:
  bool inside() const final
  {
    return true; // We expect the bases to be bound to the inside (see evaluate and post_bind).
  }

  /// \name Required by LocalUnaryIntersectionIntegrandInterface.
  /// \{

  std::unique_ptr<BaseUnaryType> copy_as_unary_intersection_integrand() const final
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalTestBasisType& test_basis, const XT::Common::Parameter& param = {}) const final
  {
    return local_dirichlet_data_->order(param) + local_diffusion_->order(param)
           + std::max(test_basis.order(param) - 1, 0);
  }

  using BaseUnaryType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const DomainType& point_in_reference_intersection,
                DynamicVector<F>& result,
                const XT::Common::Parameter& param = {}) const final
  {
    // Prepare sotrage, ...
    BaseUnaryType::ensure_size_and_clear_results(test_basis, result, param);
    // evaluate ...
    const auto point_in_inside_reference_element =
        BaseUnaryType::intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = BaseUnaryType::intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis.jacobians(point_in_inside_reference_element, test_basis_grads_, param);
    // ... data functions, ...
    const auto diffusion = local_diffusion_->evaluate(point_in_inside_reference_element, param);
    const auto dirichlet_data = local_dirichlet_data_->evaluate(point_in_inside_reference_element, param);
    // ... and finally compute the integrand.
    const size_t size = test_basis.size(param);
    for (size_t ii = 0; ii < size; ++ii)
      result[ii] += -1.0 * symmetry_prefactor_ * dirichlet_data * ((diffusion * test_basis_grads_[ii][0]) * normal);
  } // ... evaluate(...)

  /// \}
  /// \name Required by LocalBinaryIntersectionIntegrandInterface.
  /// \{

  std::unique_ptr<BaseBinaryType> copy_as_binary_intersection_integrand() const final
  {
    return std::make_unique<ThisType>(*this);
  }

  int order(const LocalTestBasisType& test_basis,
            const LocalAnsatzBasisType& ansatz_basis,
            const XT::Common::Parameter& param = {}) const final
  {
    return local_diffusion_->order(param) + test_basis.order(param) + ansatz_basis.order(param);
  }

  using BaseBinaryType::evaluate;

  void evaluate(const LocalTestBasisType& test_basis,
                const LocalAnsatzBasisType& ansatz_basis,
                const DomainType& point_in_reference_intersection,
                DynamicMatrix<F>& result,
                const XT::Common::Parameter& param = {}) const final
  {
    // Prepare sotrage, ...
    BaseBinaryType::ensure_size_and_clear_results(test_basis, ansatz_basis, result, param);
    // evaluate ...
    const auto point_in_inside_reference_element =
        BaseBinaryType::intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = BaseBinaryType::intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    test_basis.jacobians(point_in_inside_reference_element, test_basis_grads_, param);
    ansatz_basis.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    ansatz_basis.jacobians(point_in_inside_reference_element, ansatz_basis_grads_, param);
    // ... data functions, ...
    const auto diffusion = local_diffusion_->evaluate(point_in_inside_reference_element, param);
    // ... and finally compute the integrand.
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj) {
        result[ii][jj] += -1.0 * ((diffusion * ansatz_basis_grads_[jj][0]) * normal) * test_basis_values_[ii];
        result[ii][jj] +=
            -1.0 * symmetry_prefactor_ * ansatz_basis_values_[jj] * ((diffusion * test_basis_grads_[ii][0]) * normal);
      }
  } // ... evaluate(...)

  /// \}

private:
  const double symmetry_prefactor_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, d, d>> diffusion_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E>> dirichlet_data_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_diffusion_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E>::LocalFunctionType> local_dirichlet_data_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // class DirichletCoupling


} // namespace LocalLaplaceIPDGIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_LAPLACE_IPDG_HH
