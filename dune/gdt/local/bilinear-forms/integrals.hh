// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_LOCAL_BILINEAR_FORMS_INTEGRALS_HH
#define DUNE_GDT_LOCAL_BILINEAR_FORMS_INTEGRALS_HH

#include <dune/geometry/quadraturerules.hh>

#include <dune/gdt/local/integrands/interfaces.hh>
#include <dune/gdt/local/integrands/generic.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * For an explanation of the template arguments \sa LocalElementBilinearFormInterface
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalElementIntegralBilinearForm : public LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = LocalElementIntegralBilinearForm;
  using BaseType = LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using IntegrandType = LocalBinaryElementIntegrandInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;
  using GenericIntegrand = GenericLocalBinaryElementIntegrand<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

  LocalElementIntegralBilinearForm(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy())
    , over_integrate_(over_integrate)
  {}

  LocalElementIntegralBilinearForm(typename GenericIntegrand::GenericOrderFunctionType order_function,
                                   typename GenericIntegrand::GenericEvaluateFunctionType evaluate_function,
                                   const XT::Common::ParameterType& param_type = {},
                                   const int over_integrate = 0)
    : BaseType(param_type)
    , integrand_(GenericIntegrand(order_function, evaluate_function).copy())
    , over_integrate_(over_integrate)
  {}

  LocalElementIntegralBilinearForm(const ThisType& other)
    : BaseType(other.parameter_type())
    , integrand_(other.integrand_->copy())
    , over_integrate_(other.over_integrate_)
  {}

  LocalElementIntegralBilinearForm(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const LocalTestBasisType& test_basis,
              const LocalAnsatzBasisType& ansatz_basis,
              DynamicMatrix<F>& result,
              const XT::Common::Parameter& param = {}) const override final
  {
    // prepare integand
    const auto& element = ansatz_basis.element();
    assert(test_basis.element() == element && "This must not happen!");
    integrand_->bind(element);
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // loop over all quadrature points
    const auto integrand_order = integrand_->order(test_basis, ansatz_basis) + over_integrate_;
    for (const auto& quadrature_point : QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto point_in_reference_element = quadrature_point.position();
      // integration factors
      const auto integration_factor = element.geometry().integrationElement(point_in_reference_element);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_->evaluate(test_basis, ansatz_basis, point_in_reference_element, integrand_values_, param);
      assert(integrand_values_.rows() >= rows && "This must not happen!");
      assert(integrand_values_.cols() >= cols && "This must not happen!");
      // compute integral
      for (size_t ii = 0; ii < rows; ++ii)
        for (size_t jj = 0; jj < cols; ++jj)
          result[ii][jj] += integrand_values_[ii][jj] * integration_factor * quadrature_weight;
    } // loop over all quadrature points
  } // ... apply(...)

private:
  mutable std::unique_ptr<IntegrandType> integrand_;
  const int over_integrate_;
  mutable DynamicMatrix<F> integrand_values_;
}; // class LocalElementIntegralBilinearForm


/**
 * For an explanation of the template arguments \sa LocalIntersectionBilinearFormInterface
 */
template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalIntersectionIntegralBilinearForm
  : public LocalIntersectionBilinearFormInterface<I, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = LocalIntersectionIntegralBilinearForm;
  using BaseType = LocalIntersectionBilinearFormInterface<I, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using IntegrandType = LocalQuaternaryIntersectionIntegrandInterface<I, t_r, t_rC, TR, F, a_r, a_rC, AR>;

  LocalIntersectionIntegralBilinearForm(const IntegrandType& integrand, const int over_integrate = 0)
    : BaseType(integrand.parameter_type())
    , integrand_(integrand.copy())
    , over_integrate_(over_integrate)
  {}

  LocalIntersectionIntegralBilinearForm(const ThisType& other)
    : BaseType(other.parameter_type())
    , integrand_(other.integrand_->copy())
    , over_integrate_(other.over_integrate_)
  {}

  LocalIntersectionIntegralBilinearForm(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const IntersectionType& intersection,
              const LocalTestBasisType& test_basis_inside,
              const LocalAnsatzBasisType& ansatz_basis_inside,
              const LocalTestBasisType& test_basis_outside,
              const LocalAnsatzBasisType& ansatz_basis_outside,
              DynamicMatrix<F>& result_in_in,
              DynamicMatrix<F>& result_in_out,
              DynamicMatrix<F>& result_out_in,
              DynamicMatrix<F>& result_out_out,
              const XT::Common::Parameter& param = {}) const override final
  {
    // prepare integand
    integrand_->bind(intersection);
    // prepare storage
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
    // loop over all quadrature points
    const size_t integrand_order =
        integrand_->order(test_basis_inside, ansatz_basis_inside, test_basis_outside, ansatz_basis_outside)
        + over_integrate_;
    for (const auto& quadrature_point : QuadratureRules<D, d - 1>::rule(
             intersection.geometry().type(), XT::Common::numeric_cast<int>(integrand_order))) {
      const auto point_in_reference_intersection = quadrature_point.position();
      // integration factors
      const auto integration_factor = intersection.geometry().integrationElement(point_in_reference_intersection);
      const auto quadrature_weight = quadrature_point.weight();
      // evaluate the integrand
      integrand_->evaluate(test_basis_inside,
                           ansatz_basis_inside,
                           test_basis_outside,
                           ansatz_basis_outside,
                           point_in_reference_intersection,
                           integrand_values_in_in_,
                           integrand_values_in_out_,
                           integrand_values_out_in_,
                           integrand_values_out_out_,
                           param);
      assert(integrand_values_in_in_.rows() >= rows_in && "This must not happen!");
      assert(integrand_values_in_in_.cols() >= cols_in && "This must not happen!");
      assert(integrand_values_in_out_.rows() >= rows_in && "This must not happen!");
      assert(integrand_values_in_out_.cols() >= cols_out && "This must not happen!");
      assert(integrand_values_out_in_.rows() >= rows_out && "This must not happen!");
      assert(integrand_values_out_in_.cols() >= cols_in && "This must not happen!");
      assert(integrand_values_out_out_.rows() >= rows_out && "This must not happen!");
      assert(integrand_values_out_out_.cols() >= cols_out && "This must not happen!");
      // compute integral
      for (size_t ii = 0; ii < rows_in; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_in_in[ii][jj] += integrand_values_in_in_[ii][jj] * integration_factor * quadrature_weight;
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_in_out[ii][jj] += integrand_values_in_out_[ii][jj] * integration_factor * quadrature_weight;
      }
      for (size_t ii = 0; ii < rows_out; ++ii) {
        for (size_t jj = 0; jj < cols_in; ++jj)
          result_out_in[ii][jj] += integrand_values_out_in_[ii][jj] * integration_factor * quadrature_weight;
        for (size_t jj = 0; jj < cols_out; ++jj)
          result_out_out[ii][jj] += integrand_values_out_out_[ii][jj] * integration_factor * quadrature_weight;
      }
    } // loop over all quadrature points
  } // ... apply2(...)

private:
  mutable std::unique_ptr<IntegrandType> integrand_;
  const int over_integrate_;
  mutable DynamicMatrix<F> integrand_values_in_in_;
  mutable DynamicMatrix<F> integrand_values_in_out_;
  mutable DynamicMatrix<F> integrand_values_out_in_;
  mutable DynamicMatrix<F> integrand_values_out_out_;
}; // class LocalIntersectionIntegralBilinearForm


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_BILINEAR_FORMS_INTEGRALS_HH
