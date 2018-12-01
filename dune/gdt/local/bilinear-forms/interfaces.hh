// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_BILINEAR_FORMS_INTERFACES_HH
#define DUNE_GDT_LOCAL_BILINEAR_FORMS_INTERFACES_HH

#include <memory>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

namespace Dune {
namespace GDT {


/**
 * Interface for local bilinear forms associated with grid elements.
 *
 * \note Regarding SMP: the bilinear form is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Element,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class TestRangeField = double,
          class Field = double,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols,
          class AnsatzRangeField = TestRangeField>
class LocalElementBilinearFormInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalElementBilinearFormInterface<Element,
                                                     test_range_dim,
                                                     test_range_dim_cols,
                                                     TestRangeField,
                                                     Field,
                                                     ansatz_range_dim,
                                                     ansatz_range_dim_cols,
                                                     AnsatzRangeField>;

public:
  using E = Element;
  using D = typename Element::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using TR = TestRangeField;
  static const constexpr size_t t_r = test_range_dim;
  static const constexpr size_t t_rC = test_range_dim_cols;

  using AR = AnsatzRangeField;
  static const constexpr size_t a_r = ansatz_range_dim;
  static const constexpr size_t a_rC = ansatz_range_dim_cols;

  using ElementType = Element;
  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::ElementFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalElementBilinearFormInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalElementBilinearFormInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Computes the application of this bilinear form for all combinations of functions from the bases.
   */
  virtual void apply2(const LocalTestBasisType& test_basis,
                      const LocalAnsatzBasisType& ansatz_basis,
                      DynamicMatrix<F>& result,
                      const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   */
  virtual DynamicMatrix<F> apply2(const LocalTestBasisType& test_basis,
                                  const LocalAnsatzBasisType& ansatz_basis,
                                  const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result(test_basis.size(param), ansatz_basis.size(param), 0);
    this->apply2(test_basis, ansatz_basis, result, param);
    return result;
  }
}; // class LocalElementBilinearFormInterface


/**
 * Interface for local bilinear forms associated with grid intersections.
 *
 * \note Regarding SMP: the bilinear form is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Intersection,
          size_t test_range_dim = 1,
          size_t test_range_dim_cols = 1,
          class TestRangeField = double,
          class Field = double,
          size_t ansatz_range_dim = test_range_dim,
          size_t ansatz_range_dim_cols = test_range_dim_cols,
          class AnsatzRangeField = TestRangeField>
class LocalIntersectionBilinearFormInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalIntersectionBilinearFormInterface<Intersection,
                                                          test_range_dim,
                                                          test_range_dim_cols,
                                                          TestRangeField,
                                                          Field,
                                                          ansatz_range_dim,
                                                          ansatz_range_dim_cols,
                                                          AnsatzRangeField>;

public:
  using IntersectionType = Intersection;
  using ElementType = XT::Grid::extract_inside_element_t<Intersection>;

  using I = Intersection;
  using E = ElementType;
  using D = typename ElementType::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using TR = TestRangeField;
  static const constexpr size_t t_r = test_range_dim;
  static const constexpr size_t t_rC = test_range_dim_cols;

  using AR = AnsatzRangeField;
  static const constexpr size_t a_r = ansatz_range_dim;
  static const constexpr size_t a_rC = ansatz_range_dim_cols;

  using LocalTestBasisType = XT::Functions::ElementFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::ElementFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalIntersectionBilinearFormInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalIntersectionBilinearFormInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Computes the application of this bilinear form for all combinations of functions from the bases.
   */
  virtual void apply2(const IntersectionType& intersection,
                      const LocalTestBasisType& test_basis_inside,
                      const LocalAnsatzBasisType& ansatz_basis_inside,
                      const LocalTestBasisType& test_basis_outside,
                      const LocalAnsatzBasisType& ansatz_basis_outside,
                      DynamicMatrix<F>& result_in_in,
                      DynamicMatrix<F>& result_in_out,
                      DynamicMatrix<F>& result_out_in,
                      DynamicMatrix<F>& result_out_out,
                      const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   */
  virtual std::array<DynamicMatrix<F>, 4> apply2(const IntersectionType& intersection,
                                                 const LocalTestBasisType& test_basis_inside,
                                                 const LocalAnsatzBasisType& ansatz_basis_inside,
                                                 const LocalTestBasisType& test_basis_outside,
                                                 const LocalAnsatzBasisType& ansatz_basis_outside,
                                                 const XT::Common::Parameter& param = {}) const
  {
    DynamicMatrix<F> result_in_in(test_basis_inside.size(param), ansatz_basis_inside.size(param), 0);
    DynamicMatrix<F> result_in_out(test_basis_inside.size(param), ansatz_basis_outside.size(param), 0);
    DynamicMatrix<F> result_out_in(test_basis_outside.size(param), ansatz_basis_inside.size(param), 0);
    DynamicMatrix<F> result_out_out(test_basis_outside.size(param), ansatz_basis_outside.size(param), 0);
    this->apply2(intersection,
                 test_basis_inside,
                 ansatz_basis_inside,
                 test_basis_outside,
                 ansatz_basis_outside,
                 result_in_in,
                 result_in_out,
                 result_out_in,
                 result_out_out,
                 param);
    return {result_in_in, result_in_out, result_out_in, result_out_out};
  } // ... apply(...)
}; // class LocalIntersectionBilinearFormInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_BILINEAR_FORMS_INTERFACES_HH
