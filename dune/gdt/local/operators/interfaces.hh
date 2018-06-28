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

#ifndef DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
#define DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH

#include <memory>
#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template <class SourceVector,
          class SourceGridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
class LocalElementOperatorInterface : public XT::Common::ParametricInterface
{
  static_assert(
      std::is_same<XT::Grid::extract_entity_t<SourceGridView>, XT::Grid::extract_entity_t<RangeGridView>>::value, "");

public:
  using SV = SourceVector;
  using SGV = SourceGridView;
  static const constexpr size_t s_r = source_range_dim;
  static const constexpr size_t s_rC = source_range_dim_cols;
  using SR = SourceField;
  using SourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SR>;

  using RV = RangeVector;
  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_range_dim;
  static const constexpr size_t r_rC = range_range_dim_cols;
  using RR = RangeField;
  using LocalRangeType = LocalDiscreteFunction<RV, RGV, r_r, r_rC, RR>;

  using ThisType = LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SR, r_r, r_rC, RR, RGV, RV>;

  LocalElementOperatorInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalElementOperatorInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  virtual void
  apply(const SourceType& source, LocalRangeType& local_range, const XT::Common::Parameter& param = {}) const = 0;
}; // class LocalElementOperatorInterface


template <class Intersection,
          class SourceVector,
          class SourceGridView,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class InsideRangeGridView = SourceGridView,
          class InsideRangeVector = SourceVector,
          class OutsideRangeGridView = InsideRangeGridView,
          class OutsideRangeVector = InsideRangeVector>
class LocalIntersectionOperatorInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");
  static_assert(std::is_same<typename Intersection::Entity, XT::Grid::extract_entity_t<InsideRangeGridView>>::value,
                "");
  static_assert(std::is_same<typename Intersection::Entity, XT::Grid::extract_entity_t<OutsideRangeGridView>>::value,
                "");

public:
  static const constexpr size_t d = Intersection::dimension;
  using I = Intersection;
  using IntersectionType = Intersection;

  using SV = SourceVector;
  using SGV = SourceGridView;
  static const constexpr size_t s_r = source_range_dim;
  static const constexpr size_t s_rC = source_range_dim_cols;
  using SF = SourceField;
  using SourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>;

  using IRV = InsideRangeVector;
  using IRGV = InsideRangeGridView;
  static const constexpr size_t r_r = range_range_dim;
  static const constexpr size_t r_rC = range_range_dim_cols;
  using RF = RangeField;
  using LocalInsideRangeType = LocalDiscreteFunction<IRV, IRGV, r_r, r_rC, RF>;

  using ORV = OutsideRangeVector;
  using ORGV = OutsideRangeGridView;
  using LocalOutsideRangeType = LocalDiscreteFunction<ORV, ORGV, r_r, r_rC, RF>;

  using ThisType = LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, IRGV, IRV, ORGV, ORV>;

  LocalIntersectionOperatorInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalIntersectionOperatorInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
    * \note Presumes that local_range_inside is already bound to intersection.inside() and local_range_outside is
    *       already bound to intersection.outside()!
    **/
  virtual void apply(const SourceType& source,
                     const IntersectionType& intersection,
                     LocalInsideRangeType& local_range_inside,
                     LocalOutsideRangeType& local_range_outside,
                     const XT::Common::Parameter& param = {}) const = 0;
}; // class LocalIntersectionOperatorInterface


/**
 * Interface for local two-forms associated with grid elements.
 *
 * \note Regarding SMP: the two-form is copied for each thread, so
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
class LocalElementTwoFormInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

  using ThisType = LocalElementTwoFormInterface<Element,
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
  using LocalTestBasisType = XT::Functions::LocalFunctionSetInterface<E, t_r, t_rC, TR>;
  using LocalAnsatzBasisType = XT::Functions::LocalFunctionSetInterface<E, a_r, a_rC, AR>;

  LocalElementTwoFormInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~LocalElementTwoFormInterface() = default;

  virtual std::unique_ptr<ThisType> copy() const = 0;

  /**
   * Computes the application of this two-form for all combinations of functions from the bases.
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
}; // class LocalElementTwoFormInterface


#if 0
template <class TestBaseEntity,
          class Intersection,
          class AnsatzBaseEntity = TestBaseEntity,
          class TestBaseNeighbor = TestBaseEntity,
          class AnsatzBaseNeighbor = AnsatzBaseEntity,
          class Field = typename TestBaseEntity::RangeFieldType>
class LocalCouplingTwoFormInterface
{
  static_assert(XT::Functions::is_localfunction_set<TestBaseEntity>::value, "");
  static_assert(XT::Functions::is_localfunction_set<AnsatzBaseEntity>::value, "");
  static_assert(XT::Functions::is_localfunction_set<TestBaseNeighbor>::value, "");
  static_assert(XT::Functions::is_localfunction_set<AnsatzBaseNeighbor>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

public:
  typedef TestBaseEntity TestBaseEntityType;
  typedef AnsatzBaseEntity AnsatzBaseEntityType;
  typedef TestBaseNeighbor TestBaseNeighborType;
  typedef AnsatzBaseNeighbor AnsatzBaseNeighborType;
  typedef Intersection IntersectionType;
  typedef Field FieldType;

  virtual ~LocalCouplingTwoFormInterface() = default;

  virtual void apply2(const TestBaseEntityType& test_base_en,
                      const AnsatzBaseEntityType& ansatz_base_en,
                      const TestBaseNeighborType& test_base_ne,
                      const AnsatzBaseNeighborType& ansatz_base_ne,
                      const IntersectionType& intersection,
                      DynamicMatrix<FieldType>& ret_en_en,
                      DynamicMatrix<FieldType>& ret_ne_ne,
                      DynamicMatrix<FieldType>& ret_en_ne,
                      DynamicMatrix<FieldType>& ret_ne_en) const = 0;
}; // class LocalCouplingTwoFormInterface


template <class TestBase,
          class Intersection,
          class AnsatzBase = TestBase,
          class Field = typename TestBase::RangeFieldType>
class LocalBoundaryTwoFormInterface
{
  static_assert(XT::Functions::is_localfunction_set<TestBase>::value, "");
  static_assert(XT::Functions::is_localfunction_set<AnsatzBase>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

public:
  typedef TestBase TestBaseType;
  typedef AnsatzBase AnsatzBaseType;
  typedef Intersection IntersectionType;
  typedef Field FieldType;

  virtual ~LocalBoundaryTwoFormInterface() = default;

  virtual void apply2(const TestBaseType& test_base,
                      const AnsatzBaseType& ansatz_base,
                      const IntersectionType& intersection,
                      Dune::DynamicMatrix<FieldType>& ret) const = 0;

  DynamicMatrix<FieldType> apply2(const TestBaseType& test_base,
                                  const AnsatzBaseType& ansatz_base,
                                  const IntersectionType& /*intersection*/) const
  {
    DynamicMatrix<FieldType> ret(test_base.size(), ansatz_base.size(), 0.);
    apply2(test_base, ansatz_base, ret);
    return ret;
  }
}; // class LocalBoundaryTwoFormInterface
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
