// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   René Fritze     (2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014)

#warning This header is deprecated, use and include <dune/gdt/local/integrands/laplace-ipdg.hh> instead!

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH
#  define DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH

#  include <dune/xt/common/deprecated.hh>
#  include <dune/xt/functions/grid-function.hh>
#  include <dune/xt/functions/interfaces/grid-function.hh>

#  include "interfaces.hh"
#  include "ipdg.hh"

namespace Dune {
namespace GDT {

/**
 * \brief To replace the ones in LocalEllipticIpdgIntegrands at some point.
 */
namespace LocalEllipticIPDGIntegrands {


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
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I>;
  using ThisType = InnerCoupling;

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
                XT::Functions::GridFunction<E, d, d> weight_function = {1.})
    : BaseType(diffusion.parameter_type() + weight_function.parameter_type())
    , symmetry_prefactor_(symmetry_prefactor)
    , diffusion_(diffusion)
    , weight_(weight_function)
    , local_diffusion_in_(diffusion_.local_function())
    , local_diffusion_out_(diffusion_.local_function())
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  InnerCoupling(const ThisType& other)
    : BaseType(other.parameter_type())
    , symmetry_prefactor_(other.symmetry_prefactor_)
    , diffusion_(other.diffusion_)
    , weight_(other.weight_)
    , local_diffusion_in_(diffusion_.local_function())
    , local_diffusion_out_(diffusion_.local_function())
    , local_weight_in_(weight_.local_function())
    , local_weight_out_(weight_.local_function())
  {}

  InnerCoupling(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    DUNE_THROW_IF(!intersection.neighbor(),
                  Exceptions::integrand_error,
                  "This integrand cannot be used on a boundary intersection!");
    const auto inside_element = intersection.inside();
    const auto outside_element = intersection.outside();
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
            const XT::Common::Parameter& param = {}) const override final
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
  XT::Functions::GridFunction<E, d, d> diffusion_;
  XT::Functions::GridFunction<E, d, d> weight_;
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
class DirichletCoupling : public LocalQuaternaryIntersectionIntegrandInterface<I>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I>;
  using ThisType = DirichletCoupling;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  DirichletCoupling(const double& symmetry_prefactor, XT::Functions::GridFunction<E, d, d> diffusion)
    : BaseType(diffusion.parameter_type())
    , symmetry_prefactor_(symmetry_prefactor)
    , diffusion_(diffusion)
    , local_diffusion_(diffusion_.local_function())
  {}

  DirichletCoupling(const ThisType& other)
    : BaseType(other.parameter_type())
    , symmetry_prefactor_(other.symmetry_prefactor_)
    , diffusion_(other.diffusion_)
    , local_diffusion_(diffusion_.local_function())
  {}

  DirichletCoupling(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    const auto inside_element = intersection.inside();
    local_diffusion_->bind(inside_element);
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& /*test_basis_outside*/,
            const LocalAnsatzBasisType& /*ansatz_basis_outside*/,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_->order(param) + test_basis_inside.order(param) + ansatz_basis_inside.order(param);
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
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    test_basis_inside.jacobians(point_in_inside_reference_element, test_basis_grads_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    ansatz_basis_inside.jacobians(point_in_inside_reference_element, ansatz_basis_grads_, param);
    // ... data functions, ...
    const auto diffusion = local_diffusion_->evaluate(point_in_inside_reference_element, param);
    // ... and finally compute the integrand.
    const size_t rows = test_basis_inside.size(param);
    const size_t cols = ansatz_basis_inside.size(param);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj) {
        result_in_in[ii][jj] += -1.0 * ((diffusion * ansatz_basis_grads_[jj][0]) * normal) * test_basis_values_[ii];
        result_in_in[ii][jj] +=
            -1.0 * symmetry_prefactor_ * ansatz_basis_values_[jj] * ((diffusion * test_basis_grads_[ii][0]) * normal);
      }
  } // ... evaluate(...)

private:
  const double symmetry_prefactor_;
  XT::Functions::GridFunction<E, d, d> diffusion_;
  std::unique_ptr<typename XT::Functions::GridFunctionInterface<E, d, d>::LocalFunctionType> local_diffusion_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // DirichletCoupling


} // namespace LocalEllipticIPDGIntegrands


/**
 *  \brief      Contains local integrands for the family of interior penalty discontinuous Galerkin (IPDG)
 *              discretization schemes.
 *
 *              For the choice of penalization and the role of the user input see Epshteyn, Riviere (2007):
 *              "Estimation of penalty parameters for symmetric interior penalty Galerkin methods"
 *              For the choice of the weighting see Ern, Stephansen, Zunino (2007): "A discontinuous Galerkin method
 *              with weighted averages for advection-diffusion equations with locally small and anisotropic diffusivity"
 */
namespace LocalEllipticIpdgIntegrands {


enum class /*DXT_DEPRECATED_MSG("Use the LocalLaplaceIPDGIntegrands instead (10.08.2019)!")*/ Method
{
  ipdg,
  nipdg,
  sipdg,
  swipdg,
  swipdg_affine_factor,
  swipdg_affine_tensor
};


namespace internal {


template <Method method>
struct method_dependent_typename
{
  typedef void type;
};


} // namespace internal


template <Method method>
struct method_name
{
  static_assert(AlwaysFalse<typename internal::method_dependent_typename<method>::type>::value,
                "Please add a specialization for this method!");

  static std::string value()
  {
    return "";
  }
};

template <>
struct method_name<Method::ipdg>
{
  static std::string value()
  {
    return "ipdg";
  }
};

template <>
struct method_name<Method::nipdg>
{
  static std::string value()
  {
    return "nipdg";
  }
};

template <>
struct method_name<Method::sipdg>
{
  static std::string value()
  {
    return "sipdg";
  }
};

template <>
struct method_name<Method::swipdg>
{
  static std::string value()
  {
    return "swipdg";
  }
};

template <>
struct method_name<Method::swipdg_affine_factor>
{
  static std::string value()
  {
    return "swipdg_affine_factor";
  }
};

template <>
struct method_name<Method::swipdg_affine_tensor>
{
  static std::string value()
  {
    return "swipdg_affine_tensor";
  }
};


static constexpr Method default_method = Method::swipdg;


namespace internal {


/**
 * \note see Epshteyn, Riviere, 2007
 */
// DXT_DEPRECATED_MSG("Use the LocalLaplaceIPDGIntegrands instead (10.08.2019)!")
static inline double default_beta(const size_t d)
{
  return 1.0 / (d - 1.0);
}


/**
 * \note see Epshteyn, Riviere, 2007
 */
// DXT_DEPRECATED_MSG("Use the LocalLaplaceIPDGIntegrands instead (10.08.2019)!")
static inline double inner_sigma(const size_t pol_order)
{
  double sigma = 1.0;
  if (pol_order <= 1)
    sigma *= 8.0;
  else if (pol_order <= 2)
    sigma *= 20.0;
  else if (pol_order <= 3)
    sigma *= 38.0;
  else {
#  ifndef NDEBUG
#    ifndef DUNE_GDT_DISABLE_WARNINGS
    Dune::XT::Common::TimedLogger().get("gdt.local.integrands.elliptic-ipdg.inner").warn()
        << "a polynomial order of " << pol_order << " is untested!\n"
        << "  #define DUNE_GDT_DISABLE_WARNINGS to statically disable this warning\n"
        << "  or dynamically disable warnings of the TimedLogger() instance!" << std::endl;
#    endif
#  endif
    sigma *= 50.0;
  }
  return sigma;
} // ... inner_sigma(...)


/**
 * \note see Epshteyn, Riviere, 2007
 */
// DXT_DEPRECATED_MSG("Use the LocalLaplaceIPDGIntegrands instead (10.08.2019)!")
static inline double boundary_sigma(const size_t pol_order)
{
  double sigma = 1.0;
  if (pol_order <= 1)
    sigma *= 14.0;
  else if (pol_order <= 2)
    sigma *= 38.0;
  else if (pol_order <= 3)
    sigma *= 74.0;
  else {
#  ifndef NDEBUG
#    ifndef DUNE_GDT_DISABLE_WARNINGS
    Dune::XT::Common::TimedLogger().get("gdt.local.integrands.elliptic-ipdg.boundary").warn()
        << "a polynomial order of " << pol_order << " is untested!\n"
        << "  #define DUNE_GDT_DISABLE_WARNINGS to statically disable this warning\n"
        << "  or dynamically disable warnings of the TimedLogger() instance!" << std::endl;
#    endif
#  endif
    sigma *= 100.0;
  }
  return sigma;
} // ... boundary_sigma(...)


} // namespace internal


/**
 * \sa [Epshteyn, Riviere, 2007] for the meaning of beta
 */
template <class I, class F = double, Method method = default_method>
// class DXT_DEPRECATED_MSG(
//    "Use LocalLaplaceIPDGIntegrands::InnerCoupling + LocalIPDGIntegrands::InnerPenalty} instead (10.08.2019)!")
class Inner : public LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>;
  using ThisType = Inner;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  Inner(const DiffusionFactorType& diffusion_factor,
        const DiffusionTensorType& diffusion_tensor,
        const double beta = internal::default_beta(d))
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , beta_(beta)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_in_(diffusion_factor_.local_function())
    , local_diffusion_factor_out_(diffusion_factor_.local_function())
    , local_diffusion_tensor_in_(diffusion_tensor_.local_function())
    , local_diffusion_tensor_out_(diffusion_tensor_.local_function())
  {}

  Inner(const ThisType& other)
    : BaseType(other.parameter_type())
    , beta_(other.beta_)
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_in_(diffusion_factor_.local_function())
    , local_diffusion_factor_out_(diffusion_factor_.local_function())
    , local_diffusion_tensor_in_(diffusion_tensor_.local_function())
    , local_diffusion_tensor_out_(diffusion_tensor_.local_function())
  {}

  Inner(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    DUNE_THROW_IF(!intersection.neighbor(),
                  Exceptions::integrand_error,
                  "This integrand cannot be used on a boundary intersection!");
    const auto inside_element = intersection.inside();
    const auto outside_element = intersection.outside();
    local_diffusion_factor_in_->bind(inside_element);
    local_diffusion_tensor_in_->bind(inside_element);
    local_diffusion_factor_out_->bind(outside_element);
    local_diffusion_tensor_out_->bind(outside_element);
  } // ... post_bind(...)

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(local_diffusion_factor_in_->order(param), local_diffusion_factor_out_->order(param))
           + std::max(local_diffusion_tensor_in_->order(), local_diffusion_tensor_out_->order(param))
           + std::max(test_basis_inside.order(param), test_basis_outside.order(param))
           + std::max(ansatz_basis_inside.order(param), ansatz_basis_outside.order(param));
  }

  template <Method m, class Anything = void>
  struct IPDG
  {
    static_assert(AlwaysFalse<Anything>::value, "Other methods are not implemented yet!");

    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F gamma(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& /*sigma*/,
                            const R& /*gamma*/,
                            const R& /*h*/,
                            const R& /*beta*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F weight_plus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F weight_minus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }
  }; // struct IPDG<...>

  template <class Anything>
  struct IPDG<Method::swipdg, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& diffusion_ne,
                               const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_ne * normal);
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& diffusion_en,
                                const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_en * normal);
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (sigma * gamma) / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_factor, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_ne,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_tensor_ne * normal);
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_en,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_tensor_en * normal);
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& diffusion_factor_en,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& diffusion_factor_ne,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (0.5 * (diffusion_factor_en + diffusion_factor_ne) * sigma * gamma) / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg_affine_factor, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_tensor, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& diffusion_factor_ne,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      return diffusion_factor_ne;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& diffusion_factor_en,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      return diffusion_factor_en;
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_en,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_ne,
                            const FieldVector<R, d>& normal,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (normal * (((diffusion_tensor_en + diffusion_tensor_ne) * 0.5) * normal) * sigma * gamma)
             / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg_affine_tensor, ...>

  template <class Anything>
  struct IPDG<Method::sipdg, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F gamma(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& /*gamma*/,
                            const R& h,
                            const R& beta)
    {
      return sigma / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 0.5;
    }

    template <class R>
    static inline F weight_minus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 0.5;
    }
  }; // struct IPDG<Method::sipdg, ...>

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
    // prepare sotrage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out);
    ensure_size_and_clear(result_out_in, rows_out, cols_in);
    ensure_size_and_clear(result_out_out, rows_out, cols_out);
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
    // ... data functions
    const auto diffusion_factor_in = local_diffusion_factor_in_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_tensor_in = local_diffusion_tensor_in_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_factor_out = local_diffusion_factor_out_->evaluate(point_in_outside_reference_element, param);
    const auto diffusion_tensor_out = local_diffusion_tensor_out_->evaluate(point_in_outside_reference_element, param);
    const auto diffusion_in = diffusion_tensor_in * diffusion_factor_in;
    const auto diffusion_out = diffusion_tensor_out * diffusion_factor_out;
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(test_basis_inside.order(param),
                 std::max(ansatz_basis_inside.order(param),
                          std::max(test_basis_outside.order(param), ansatz_basis_outside.order(param))));
    const F sigma = internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const F delta_plus = IPDG<method>::delta_plus(diffusion_factor_out, diffusion_tensor_out, diffusion_out, normal);
    const F delta_minus = IPDG<method>::delta_minus(diffusion_factor_in, diffusion_tensor_in, diffusion_in, normal);
    const F gamma = IPDG<method>::gamma(delta_plus, delta_minus);
    const F penalty = IPDG<method>::penalty(diffusion_factor_in,
                                            diffusion_tensor_out,
                                            diffusion_factor_out,
                                            diffusion_tensor_in,
                                            normal,
                                            sigma,
                                            gamma,
                                            this->intersection().geometry().volume(),
                                            beta_);
    const F weight_plus = IPDG<method>::weight_plus(delta_plus, delta_minus);
    const F weight_minus = IPDG<method>::weight_minus(delta_plus, delta_minus);
    // compute integrand
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj) {
        // consistency term
        result_in_in[ii][jj] +=
            -weight_minus * ((diffusion_in * ansatz_basis_in_grads_[jj][0]) * normal) * test_basis_in_values_[ii];
        // symmetry term
        result_in_in[ii][jj] +=
            -weight_minus * ansatz_basis_in_values_[jj] * ((diffusion_in * test_basis_in_grads_[ii][0]) * normal);
        // penalty term
        result_in_in[ii][jj] += penalty * ansatz_basis_in_values_[jj] * test_basis_in_values_[ii];
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        // consistency term
        result_in_out[ii][jj] +=
            -weight_plus * ((diffusion_out * ansatz_basis_out_grads_[jj][0]) * normal) * test_basis_in_values_[ii];
        // symmetry term
        result_in_out[ii][jj] +=
            weight_minus * ansatz_basis_out_values_[jj] * ((diffusion_in * test_basis_in_grads_[ii][0]) * normal);
        // penalty term
        result_in_out[ii][jj] += -1.0 * penalty * ansatz_basis_out_values_[jj] * test_basis_in_values_[ii];
      }
    }
    for (size_t ii = 0; ii < rows_out; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj) {
        // consistency term
        result_out_in[ii][jj] +=
            weight_minus * ((diffusion_in * ansatz_basis_in_grads_[jj][0]) * normal) * test_basis_out_values_[ii];
        // symmetry term
        result_out_in[ii][jj] +=
            -weight_plus * ansatz_basis_in_values_[jj] * ((diffusion_out * test_basis_out_grads_[ii][0]) * normal);
        // penalty term
        result_out_in[ii][jj] += -1.0 * penalty * ansatz_basis_in_values_[jj] * test_basis_out_values_[ii];
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        // consistency term
        result_out_out[ii][jj] +=
            weight_plus * ((diffusion_out * ansatz_basis_out_grads_[jj][0]) * normal) * test_basis_out_values_[ii];
        // symmetry term
        result_out_out[ii][jj] +=
            weight_plus * ansatz_basis_out_values_[jj] * ((diffusion_out * test_basis_out_grads_[ii][0]) * normal);
        // penalty term
        result_out_out[ii][jj] += penalty * ansatz_basis_out_values_[jj] * test_basis_out_values_[ii];
      }
    }
  } // ... evaluate(...)

private:
  const double beta_;
  const DiffusionFactorType& diffusion_factor_; // These are just required ...
  const DiffusionTensorType& diffusion_tensor_; //                         ... for the copy ctor atm.
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_in_;
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_out_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_in_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_in_grads_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_out_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_in_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_out_grads_;
}; // Inner


/**
 * \sa [Epshteyn, Riviere, 2007] for the meaning of beta
 */
template <class I, class F = double, Method method = default_method>
class /*DXT_DEPRECATED_MSG(
    "Use LocalLaplaceIPDGIntegrands::DirichletCoupling + LocalIPDGIntegrands::boundaryPenalty instead (10.08.2019)!")*/
    DirichletBoundaryLhs : public LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>;
  using ThisType = DirichletBoundaryLhs;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  DirichletBoundaryLhs(const DiffusionFactorType& diffusion_factor,
                       const DiffusionTensorType& diffusion_tensor,
                       const double beta = internal::default_beta(d))
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , beta_(beta)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {}

  DirichletBoundaryLhs(const ThisType& other)
    : BaseType(other.parameter_type())
    , beta_(other.beta_)
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {}

  DirichletBoundaryLhs(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    const auto inside_element = intersection.inside();
    local_diffusion_factor_->bind(inside_element);
    local_diffusion_tensor_->bind(inside_element);
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& /*test_basis_outside*/,
            const LocalAnsatzBasisType& /*ansatz_basis_outside*/,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_factor_->order(param) + local_diffusion_tensor_->order(param)
           + test_basis_inside.order(param) + ansatz_basis_inside.order(param);
  }

  template <Method m, class Anything = void>
  struct IPDG
  {
    static_assert(AlwaysFalse<Anything>::value, "Other methods are not implemented yet!");

    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/, const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F penalty(const R& /*sigma*/, const R& /*gamma*/, const R& /*h*/, const R& /*beta*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }
  }; // struct IPDG<...>

  template <class Anything>
  struct IPDG<Method::swipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& diffusion, const FieldVector<R, d>& normal)
    {
      return normal * (diffusion * normal);
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& gamma, const R& h, const R& beta)
    {
      return (sigma * gamma) / std::pow(h, beta);
    }
  }; // struct IPDG<Method::swipdg, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_factor, Anything> : public IPDG<Method::swipdg, Anything>
  {};

  template <class Anything>
  struct IPDG<Method::swipdg_affine_tensor, Anything> : public IPDG<Method::swipdg, Anything>
  {};

  template <class Anything>
  struct IPDG<Method::sipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/, const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& /*gamma*/, const R& h, const R& beta)
    {
      return sigma / std::pow(h, beta);
    }
  }; // struct IPDG<Method::sipdg, ...>

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
    // prepare sotrage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out); // This is on purpose ...
    ensure_size_and_clear(result_out_in, rows_out, cols_in); // ... (including the resize), otherwise ...
    ensure_size_and_clear(result_out_out, rows_out, cols_out); // ... we could not use this in the standard assembler.
    // evaluate ...
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    test_basis_inside.jacobians(point_in_inside_reference_element, test_basis_grads_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    ansatz_basis_inside.jacobians(point_in_inside_reference_element, ansatz_basis_grads_, param);
    // ... data functions
    const auto diffusion_factor = local_diffusion_factor_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_tensor = local_diffusion_tensor_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion = diffusion_tensor * diffusion_factor;
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(test_basis_inside.order(param), ansatz_basis_inside.order(param));
    const F sigma = internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const F gamma = IPDG<method>::gamma(diffusion, normal);
    const F penalty = IPDG<method>::penalty(sigma, gamma, this->intersection().geometry().volume(), beta_);
    // compute integrand
    for (size_t ii = 0; ii < rows_in; ++ii) {
      for (size_t jj = 0; jj < cols_in; ++jj) {
        // consistency term
        result_in_in[ii][jj] += -1.0 * ((diffusion * ansatz_basis_grads_[jj][0]) * normal) * test_basis_values_[ii];
        // symmetry term
        result_in_in[ii][jj] += -1.0 * ansatz_basis_values_[jj] * ((diffusion * test_basis_grads_[ii][0]) * normal);
        // penalty term
        result_in_in[ii][jj] += penalty * ansatz_basis_values_[jj] * test_basis_values_[ii];
      }
    }
  } // ... evaluate(...)

private:
  const double beta_;
  const DiffusionFactorType& diffusion_factor_; // These are just required ...
  const DiffusionTensorType& diffusion_tensor_; //                         ... for the copy ctor atm.
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // DirichletBoundaryLhs


#  if 0
template <class DirichletImp, class DiffusionFactorImp, class DiffusionTensorImp, Method method>
class BoundaryRHS : public LocalFaceIntegrandInterface<internal::BoundaryRHSTraits<DirichletImp,
                                                                                   DiffusionFactorImp,
                                                                                   DiffusionTensorImp,
                                                                                   method>,
                                                       1>
{
  typedef LocalEllipticIntegrand<DiffusionFactorImp, DiffusionTensorImp> EllipticType;
  typedef BoundaryRHS<DirichletImp, DiffusionFactorImp, DiffusionTensorImp, method> ThisType;

public:
  typedef internal::BoundaryRHSTraits<DirichletImp, DiffusionFactorImp, DiffusionTensorImp, method> Traits;
  typedef typename Traits::DirichletType DirichletType;
  typedef typename Traits::DiffusionFactorType DiffusionFactorType;
  typedef typename Traits::DiffusionTensorType DiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  static const size_t d = Traits::d;

  BoundaryRHS(const DirichletType& dirichlet,
              const DiffusionFactorType& diffusion_factor,
              const DiffusionTensorType& diffusion_tensor,
              const double beta = internal::default_beta(d))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion_factor, diffusion_tensor)
    , beta_(beta)
  {
  }

  BoundaryRHS(const DirichletType& dirichlet,
              const DiffusionFactorType& diffusion_factor,
              const double beta = internal::default_beta(d))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion_factor)
    , beta_(beta)
  {
  }

  template <
      typename DiffusionType // This disables the ctor if d == 1, since factor and tensor are then identical
      ,
      typename = typename std::enable_if<(std::is_same<DiffusionType, DiffusionTensorType>::value) // and the ctors
                                         && (d > 1)
                                         && sizeof(DiffusionType)>::type> // ambiguous.
  BoundaryRHS(const DirichletType& dirichlet,
              const DiffusionType& diffusion,
              const double beta = internal::default_beta(d))
    : dirichlet_(dirichlet)
    , elliptic_(diffusion)
    , beta_(beta)
  {
  }

  BoundaryRHS(const ThisType& other) = default;
  BoundaryRHS(ThisType&& source) = default;

  /// \name Required by LocalFaceIntegrandInterface< ..., 1 >.
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(dirichlet_.local_function(entity),
                           elliptic_.diffusion_factor().local_function(entity),
                           elliptic_.diffusion_tensor().local_function(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class R, size_t r, size_t rC>
  size_t order(
      const LocalfunctionTupleType& local_functions,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, d, R, r, rC>& test_base) const
  {
    const auto local_dirichlet = std::get<0>(local_functions);
    const auto local_diffusion_factor = std::get<1>(local_functions);
    const auto local_diffusion_tensor = std::get<2>(local_functions);
    return order(*local_dirichlet, *local_diffusion_factor, *local_diffusion_tensor, test_base);
  }

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class IntersectionType, class R, size_t r, size_t rC>
  void
  evaluate(const LocalfunctionTupleType& local_functions,
           const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, d, R, r, rC>& test_base,
           const IntersectionType& intersection,
           const Dune::FieldVector<DomainFieldType, d - 1>& local_point,
           Dune::DynamicVector<R>& ret) const
  {
    const auto local_dirichlet = std::get<0>(local_functions);
    const auto local_diffusion_factor = std::get<1>(local_functions);
    const auto local_diffusion_tensor = std::get<2>(local_functions);
    evaluate(
        *local_dirichlet, *local_diffusion_factor, *local_diffusion_tensor, test_base, intersection, local_point, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rDF, size_t rCDF, size_t rDT, size_t rCDT, size_t rLR, size_t rCLR, size_t rT, size_t rCT>
  size_t order(const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, rLR, rCLR>&
                   local_dirichlet,
               const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, rDF, rCDF>&
                   local_diffusion_factor,
               const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, rDT, rCDT>&
                   local_diffusion_tensor,
               const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, d, R, rT, rCT>&
                   test_base) const
  {
    const size_t test_order = test_base.order();
    const size_t test_gradient_order = std::max(ssize_t(test_order) - 1, ssize_t(0));
    const size_t diffusionOrder = local_diffusion_factor.order() + local_diffusion_tensor.order();
    const size_t dirichletOrder = local_dirichlet.order();
    return std::max(test_order + dirichletOrder, diffusionOrder + test_gradient_order + dirichletOrder);
  } // ... order(...)

private:
  template <Method m, class Anything = void>
  struct IPDG
  {
    static_assert(AlwaysFalse<Anything>::value, "Other methods are not implemented yet!");

    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/,
                          const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F penalty(const R& /*sigma*/, const R& /*gamma*/, const R& /*h*/, const R& /*beta*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }
  }; // struct IPDG<...>

  template <class Anything>
  struct IPDG<Method::swipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& diffusion,
                          const FieldVector<R, d>& normal)
    {
      return normal * (diffusion * normal);
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& gamma, const R& h, const R& beta)
    {
      return (sigma * gamma) / std::pow(h, beta);
    }
  }; // struct IPDG<Method::swipdg, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_factor, Anything> : public IPDG<Method::swipdg, Anything>
  {
  };

  template <class Anything>
  struct IPDG<Method::swipdg_affine_tensor, Anything> : public IPDG<Method::swipdg, Anything>
  {
  };

  template <class Anything>
  struct IPDG<Method::sipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/,
                          const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& /*gamma*/, const R& h, const R& beta)
    {
      return sigma / std::pow(h, beta);
    }
  }; // struct IPDG<Method::sipdg, ...>

public:
  /// \}
  /// \name Actual implementation of evaluate.
  /// \{

  template <class R, class IntersectionType>
  void evaluate(
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, 1, 1>& local_dirichlet,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, 1, 1>&
          local_diffusion_factor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, d, R, d, d>&
          local_diffusion_tensor,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, d, R, 1, 1>& test_base,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, d - 1>& local_point,
      Dune::DynamicVector<R>& ret) const
  {
    typedef XT::Common::FieldMatrix<R, d, d> TensorType;
    // clear ret
    ret *= 0.0;
    // get local point (which is in intersection coordinates) in entity coordinates
    const auto local_point_entity = intersection.geometryInInside().global(local_point);
    const auto normal = intersection.unitOuterNormal(local_point);
    // evaluate local functions
    const auto dirichlet_value = local_dirichlet.evaluate(local_point_entity);
    const auto diffusion_factor_value = local_diffusion_factor.evaluate(local_point_entity);
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point_entity);
    const auto diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t polorder = test_base.order();
    const F sigma = internal::boundary_sigma(polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const F gamma = IPDG<method>::gamma(diffusion_value, normal);
    const F penalty = IPDG<method>::penalty(sigma, gamma, intersection.geometry().volume(), beta_);
    // evaluate basis
    const size_t size = test_base.size();
    const auto test_values = test_base.evaluate(local_point_entity);
    const auto test_gradients = test_base.jacobian(local_point_entity);
    // compute
    assert(ret.size() >= size);
    // loop over all test basis functions
    for (size_t ii = 0; ii < size; ++ii) {
      // symmetry term
      ret[ii] += -1.0 * dirichlet_value * ((diffusion_value * test_gradients[ii][0]) * normal);
      // penalty term
      ret[ii] += penalty * dirichlet_value * test_values[ii];
    } // loop over all test basis functions
  } // ... evaluate(...)

  /// \}

  const DirichletType& dirichlet_;
  const EllipticType elliptic_;
  const double beta_;
}; // class BoundaryRHS
#  endif // 0

/**
 * \sa [Epshteyn, Riviere, 2007] for the meaning of beta
 */
template <class I, class F = double, Method method = default_method>
class /*DXT_DEPRECATED_MSG("Use LocalIPDGIntegrands::InnerPenalty instead (05.08.2019)!")*/ InnerOnlyPenalty
  : public LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>;
  using ThisType = InnerOnlyPenalty;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  InnerOnlyPenalty(const DiffusionFactorType& diffusion_factor,
                   const DiffusionTensorType& diffusion_tensor,
                   const double beta = internal::default_beta(d))
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , beta_(beta)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_in_(diffusion_factor_.local_function())
    , local_diffusion_factor_out_(diffusion_factor_.local_function())
    , local_diffusion_tensor_in_(diffusion_tensor_.local_function())
    , local_diffusion_tensor_out_(diffusion_tensor_.local_function())
  {}

  InnerOnlyPenalty(const ThisType& other)
    : BaseType(other.parameter_type())
    , beta_(other.beta_)
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_in_(diffusion_factor_.local_function())
    , local_diffusion_factor_out_(diffusion_factor_.local_function())
    , local_diffusion_tensor_in_(diffusion_tensor_.local_function())
    , local_diffusion_tensor_out_(diffusion_tensor_.local_function())
  {}

  InnerOnlyPenalty(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    DUNE_THROW_IF(!intersection.neighbor(),
                  Exceptions::integrand_error,
                  "This integrand cannot be used on a boundary intersection!");
    const auto inside_element = intersection.inside();
    const auto outside_element = intersection.outside();
    local_diffusion_factor_in_->bind(inside_element);
    local_diffusion_tensor_in_->bind(inside_element);
    local_diffusion_factor_out_->bind(outside_element);
    local_diffusion_tensor_out_->bind(outside_element);
  } // ... post_bind(...)

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& test_basis_outside,
            const LocalAnsatzBasisType& ansatz_basis_outside,
            const XT::Common::Parameter& param = {}) const override final
  {
    return std::max(local_diffusion_factor_in_->order(param), local_diffusion_factor_out_->order(param))
           + std::max(local_diffusion_tensor_in_->order(), local_diffusion_tensor_out_->order(param))
           + std::max(test_basis_inside.order(param), test_basis_outside.order(param))
           + std::max(ansatz_basis_inside.order(param), ansatz_basis_outside.order(param));
  }

private:
  template <Method m, class Anything = void>
  struct IPDG
  {
    static_assert(AlwaysFalse<Anything>::value, "Other methods are not implemented yet!");

    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F gamma(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& /*sigma*/,
                            const R& /*gamma*/,
                            const R& /*h*/,
                            const R& /*beta*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F weight_plus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F weight_minus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }
  }; // struct IPDG<...>

  template <class Anything>
  struct IPDG<Method::swipdg, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& diffusion_ne,
                               const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_ne * normal);
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& diffusion_en,
                                const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_en * normal);
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (sigma * gamma) / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_factor, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_ne,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_tensor_ne * normal);
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_en,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& normal)
    {
      return normal * (diffusion_tensor_en * normal);
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& diffusion_factor_en,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& diffusion_factor_ne,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (0.5 * (diffusion_factor_en + diffusion_factor_ne) * sigma * gamma) / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg_affine_factor, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_tensor, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& diffusion_factor_ne,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      return diffusion_factor_ne;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& diffusion_factor_en,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      return diffusion_factor_en;
    }

    template <class R>
    static inline F gamma(const R& delta_plus, const R& delta_minus)
    {
      return (delta_plus * delta_minus) / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_en,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& diffusion_tensor_ne,
                            const FieldVector<R, d>& normal,
                            const R& sigma,
                            const R& gamma,
                            const R& h,
                            const R& beta)
    {
      return (normal * (((diffusion_tensor_en + diffusion_tensor_ne) * 0.5) * normal) * sigma * gamma)
             / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& delta_plus, const R& delta_minus)
    {
      return delta_minus / (delta_plus + delta_minus);
    }

    template <class R>
    static inline F weight_minus(const R& delta_plus, const R& delta_minus)
    {
      return delta_plus / (delta_plus + delta_minus);
    }
  }; // struct IPDG<Method::swipdg_affine_tensor, ...>

  template <class Anything>
  struct IPDG<Method::sipdg, Anything>
  {
    template <class R>
    static inline F delta_plus(const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                               const XT::Common::FieldMatrix<R, d, d>& /*diffusion_ne*/,
                               const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F delta_minus(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                                const XT::Common::FieldMatrix<R, d, d>& /*diffusion_en*/,
                                const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F gamma(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F penalty(const FieldVector<R, 1>& /*diffusion_factor_en*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_en*/,
                            const FieldVector<R, 1>& /*diffusion_factor_ne*/,
                            const XT::Common::FieldMatrix<R, d, d>& /*diffusion_tensor_ne*/,
                            const FieldVector<R, d>& /*normal*/,
                            const R& sigma,
                            const R& /*gamma*/,
                            const R& h,
                            const R& beta)
    {
      return sigma / std::pow(h, beta);
    }

    template <class R>
    static inline F weight_plus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 0.5;
    }

    template <class R>
    static inline F weight_minus(const R& /*delta_plus*/, const R& /*delta_minus*/)
    {
      return 0.5;
    }
  }; // struct IPDG<Method::sipdg, ...>

public:
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
    // prepare sotrage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out);
    ensure_size_and_clear(result_out_in, rows_out, cols_in);
    ensure_size_and_clear(result_out_out, rows_out, cols_out);
    // evaluate ...
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    const auto point_in_outside_reference_element =
        this->intersection().geometryInOutside().global(point_in_reference_intersection);
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_in_values_, param);
    test_basis_outside.evaluate(point_in_outside_reference_element, test_basis_out_values_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_in_values_, param);
    ansatz_basis_outside.evaluate(point_in_outside_reference_element, ansatz_basis_out_values_, param);
    // ... data functions
    const auto diffusion_factor_in = local_diffusion_factor_in_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_tensor_in = local_diffusion_tensor_in_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_factor_out = local_diffusion_factor_out_->evaluate(point_in_outside_reference_element, param);
    const auto diffusion_tensor_out = local_diffusion_tensor_out_->evaluate(point_in_outside_reference_element, param);
    const auto diffusion_in = diffusion_tensor_in * diffusion_factor_in;
    const auto diffusion_out = diffusion_tensor_out * diffusion_factor_out;
    // compute penalty factor (see Epshteyn, Riviere, 2007)
    const size_t max_polorder =
        std::max(test_basis_inside.order(param),
                 std::max(ansatz_basis_inside.order(param),
                          std::max(test_basis_outside.order(param), ansatz_basis_outside.order(param))));
    const F sigma = internal::inner_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const F delta_plus = IPDG<method>::delta_plus(diffusion_factor_out, diffusion_tensor_out, diffusion_out, normal);
    const F delta_minus = IPDG<method>::delta_minus(diffusion_factor_in, diffusion_tensor_in, diffusion_in, normal);
    const F gamma = IPDG<method>::gamma(delta_plus, delta_minus);
    const F penalty = IPDG<method>::penalty(diffusion_factor_in,
                                            diffusion_tensor_out,
                                            diffusion_factor_out,
                                            diffusion_tensor_in,
                                            normal,
                                            sigma,
                                            gamma,
                                            this->intersection().geometry().volume(),
                                            beta_);
    // compute integrand
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
  const double beta_;
  const DiffusionFactorType& diffusion_factor_; // These are just required ...
  const DiffusionTensorType& diffusion_tensor_; //                         ... for the copy ctor atm.
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_in_;
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_out_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_in_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_out_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_in_values_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_out_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_in_values_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_out_values_;
}; //  InnerOnlyPenalty


/**
 * \sa [Epshteyn, Riviere, 2007] for the meaning of beta
 */
template <class I, class F = double, Method method = default_method>
class /*DXT_DEPRECATED_MSG("Use LocalIPDGIntegrands::BoundaryPenalty instead (05.08.2019)!")*/
    DirichletBoundaryLhsOnlyPenalty : public LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>
{
  using BaseType = LocalQuaternaryIntersectionIntegrandInterface<I, 1, 1, F, F, 1, 1, F>;
  using ThisType = DirichletBoundaryLhsOnlyPenalty;

public:
  using BaseType::d;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using DiffusionFactorType = XT::Functions::GridFunctionInterface<E, 1, 1, F>;
  using DiffusionTensorType = XT::Functions::GridFunctionInterface<E, d, d, F>;

  DirichletBoundaryLhsOnlyPenalty(const DiffusionFactorType& diffusion_factor,
                                  const DiffusionTensorType& diffusion_tensor,
                                  const double beta = internal::default_beta(d))
    : BaseType(diffusion_factor.parameter_type() + diffusion_tensor.parameter_type())
    , beta_(beta)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {}

  DirichletBoundaryLhsOnlyPenalty(const ThisType& other)
    : BaseType(other.parameter_type())
    , beta_(other.beta_)
    , diffusion_factor_(other.diffusion_factor_)
    , diffusion_tensor_(other.diffusion_tensor_)
    , local_diffusion_factor_(diffusion_factor_.local_function())
    , local_diffusion_tensor_(diffusion_tensor_.local_function())
  {}

  DirichletBoundaryLhsOnlyPenalty(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

protected:
  void post_bind(const IntersectionType& intersection) override final
  {
    const auto inside_element = intersection.inside();
    local_diffusion_factor_->bind(inside_element);
    local_diffusion_tensor_->bind(inside_element);
  }

public:
  int order(const LocalTestBasisType& test_basis_inside,
            const LocalAnsatzBasisType& ansatz_basis_inside,
            const LocalTestBasisType& /*test_basis_outside*/,
            const LocalAnsatzBasisType& /*ansatz_basis_outside*/,
            const XT::Common::Parameter& param = {}) const override final
  {
    return local_diffusion_factor_->order(param) + local_diffusion_tensor_->order(param)
           + test_basis_inside.order(param) + ansatz_basis_inside.order(param);
  }

private:
  template <Method m, class Anything = void>
  struct IPDG
  {
    static_assert(AlwaysFalse<Anything>::value, "Other methods are not implemented yet!");

    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/, const FieldVector<R, d>& /*normal*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }

    template <class R>
    static inline F penalty(const R& /*sigma*/, const R& /*gamma*/, const R& /*h*/, const R& /*beta*/)
    {
      static_assert(AlwaysFalse<R>::value, "Other methods are not implemented yet!");
      return 0.;
    }
  }; // struct IPDG<...>

  template <class Anything>
  struct IPDG<Method::swipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& diffusion, const FieldVector<R, d>& normal)
    {
      return normal * (diffusion * normal);
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& gamma, const R& h, const R& beta)
    {
      return (sigma * gamma) / std::pow(h, beta);
    }
  }; // struct IPDG<Method::swipdg, ...>

  template <class Anything>
  struct IPDG<Method::swipdg_affine_factor, Anything> : public IPDG<Method::swipdg, Anything>
  {};

  template <class Anything>
  struct IPDG<Method::swipdg_affine_tensor, Anything> : public IPDG<Method::swipdg, Anything>
  {};

  template <class Anything>
  struct IPDG<Method::sipdg, Anything>
  {
    template <class R>
    static inline F gamma(const XT::Common::FieldMatrix<R, d, d>& /*diffusion*/, const FieldVector<R, d>& /*normal*/)
    {
      return 1.0;
    }

    template <class R>
    static inline F penalty(const R& sigma, const R& /*gamma*/, const R& h, const R& beta)
    {
      return sigma / std::pow(h, beta);
    }
  }; // struct IPDG<Method::sipdg, ...>

public:
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
    // prepare sotrage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out); // This is on purpose ...
    ensure_size_and_clear(result_out_in, rows_out, cols_in); // ... (including the resize), otherwise ...
    ensure_size_and_clear(result_out_out, rows_out, cols_out); // ... we could not use this in the standard assembler.
    // evaluate ...
    const auto point_in_inside_reference_element =
        this->intersection().geometryInInside().global(point_in_reference_intersection);
    const auto normal = this->intersection().unitOuterNormal(point_in_reference_intersection);
    // ... basis functions and ...
    test_basis_inside.evaluate(point_in_inside_reference_element, test_basis_values_, param);
    ansatz_basis_inside.evaluate(point_in_inside_reference_element, ansatz_basis_values_, param);
    // ... data functions
    const auto diffusion_factor = local_diffusion_factor_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion_tensor = local_diffusion_tensor_->evaluate(point_in_inside_reference_element, param);
    const auto diffusion = diffusion_tensor * diffusion_factor;
    // compute penalty (see Epshteyn, Riviere, 2007)
    const size_t max_polorder = std::max(test_basis_inside.order(param), ansatz_basis_inside.order(param));
    const F sigma = internal::boundary_sigma(max_polorder);
    // compute weighting (see Ern, Stephansen, Zunino 2007)
    const F gamma = IPDG<method>::gamma(diffusion, normal);
    const F penalty = IPDG<method>::penalty(sigma, gamma, this->intersection().geometry().volume(), beta_);
    // compute integrand
    for (size_t ii = 0; ii < rows_in; ++ii)
      for (size_t jj = 0; jj < cols_in; ++jj)
        result_in_in[ii][jj] += penalty * ansatz_basis_values_[jj] * test_basis_values_[ii];
  } // ... evaluate(...)

private:
  const double beta_;
  const DiffusionFactorType& diffusion_factor_; // These are just required ...
  const DiffusionTensorType& diffusion_tensor_; //                         ... for the copy ctor atm.
  std::unique_ptr<typename DiffusionFactorType::LocalFunctionType> local_diffusion_factor_;
  std::unique_ptr<typename DiffusionTensorType::LocalFunctionType> local_diffusion_tensor_;
  mutable std::vector<typename LocalTestBasisType::RangeType> test_basis_values_;
  mutable std::vector<typename LocalTestBasisType::DerivativeRangeType> test_basis_grads_;
  mutable std::vector<typename LocalAnsatzBasisType::RangeType> ansatz_basis_values_;
  mutable std::vector<typename LocalAnsatzBasisType::DerivativeRangeType> ansatz_basis_grads_;
}; // DirichletBoundaryLhsOnlyPenalty


} // namespace LocalEllipticIpdgIntegrands
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ELLIPTIC_IPDG_HH
