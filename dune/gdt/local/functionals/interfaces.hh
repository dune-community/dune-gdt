// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/exceptions.hh>


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
}; // class LocalElementFunctionalInterface


/**
 * Interface for local functionals associated with grid intersections.
 *
 * \note Although the apply gets two arguments, this is a functional and not a bilinear form (the two arguments are not
 *       ansatz and test basis on the intersection, but the test basis on the inside and on the outisde, respectively)!
 * \note Regarding SMP: the functional is copied for each thread, so
 *       - no shared mutable state between copies to be thread safe, but
 *       - local mutable state is ok.
 */
template <class Intersection,
          size_t range_dim = 1,
          size_t range_dim_cols = 1,
          class RangeField = double,
          class Field = double>
class LocalIntersectionFunctionalInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

  using ThisType = LocalIntersectionFunctionalInterface<Intersection, range_dim, range_dim_cols, RangeField>;

public:
  using IntersectionType = Intersection;
  using ElementType = typename Intersection::Entity;
  using D = typename ElementType::Geometry::ctype;
  static const constexpr size_t d = ElementType::dimension;

  using I = IntersectionType;
  using E = ElementType;
  using F = Field;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_cols;

  using LocalBasisType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

  LocalIntersectionFunctionalInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {}

  virtual ~LocalIntersectionFunctionalInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Flag to document which element the basis is expected to be bound to.
   */
  virtual bool inside() const
  {
    return true;
  }

  /**
   * Computes the application of this functional for each basis function.
   */
  virtual void apply(const IntersectionType& intersection,
                     const LocalBasisType& test_basis,
                     DynamicVector<F>& result,
                     const XT::Common::Parameter& param = {}) const = 0;

  /**
   * This method is provided for convenience and should not be used within library code.
   */
  virtual DynamicVector<F> apply(const IntersectionType& intersection,
                                 const LocalBasisType& test_basis,
                                 const XT::Common::Parameter& param = {}) const
  {
    DynamicVector<F> ret(test_basis.size(param), 0);
    this->apply(intersection, test_basis, ret, param);
    return ret;
  }
}; // class LocalIntersectionFunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
