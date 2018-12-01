// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>


namespace Dune {
namespace GDT {


/**
 * Interface for local functionals associated with grid elements.
 *
 * \note Regarding SMP: the functional is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Element,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double,
          class Field = double>
class LocalElementFunctionalInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalElementFunctionalInterface<Element, range_dim, range_dim_cols, RangeField>;

public:
  using E = Element;
  using D = typename Element::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using F = Field;

  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using ElementType = Element;
  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

  LocalElementFunctionalInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalElementFunctionalInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Computes the application of this functional to each function of the basis.
   */
  virtual void
  apply(const LocalBasisType& basis, DynamicVector<F>& result, const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   */
  virtual DynamicVector<F> apply(const LocalBasisType& basis, const XT::Common::Parameter& param = {}) const
  {
    DynamicVector<F> ret(basis.size(param));
    this->apply(basis, ret, param);
    return ret;
  }
}; // class LocalFunctionalInterface


#if 0
template <class TestBase, class Intersection, class Field = typename TestBase::RangeFieldType>
class LocalFaceFunctionalInterface
{
  static_assert(XT::Functions::is_localfunction_set<TestBase>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

public:
  typedef TestBase TestBaseType;
  typedef Intersection IntersectionType;
  typedef Field FieldType;

  virtual ~LocalFaceFunctionalInterface() = default;

  virtual void
  apply(const TestBaseType& test_basis, const IntersectionType& intersection, DynamicVector<FieldType>& ret) const = 0;

  DynamicVector<FieldType> apply(const TestBaseType& test_basis, const IntersectionType& intersection) const
  {
    DynamicVector<FieldType> ret(test_basis.size(), 0.);
    apply(test_basis, intersection, ret);
    return ret;
  }
}; // class LocalFaceFunctionalInterface
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
