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

#include <vector>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/local/discretefunction.hh>

namespace Dune {
namespace GDT {
#if 0
namespace internal {

/// \todo drop and implement the is_... properly
class IsLocalOperator
{
};

/// \todo drop and implement the is_... properly
class IsLocalCouplingOperator
{
};

/// \todo drop and implement the is_... properly
class IsLocalBoundaryOperator
{
};


} // namespace internal


template <class Traits>
class LocalOperatorInterface : public XT::CRTPInterface<LocalOperatorInterface<Traits>, Traits>,
                               internal::IsLocalOperator
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   * \brief Applies the local operator.
   * \param source Should be a XT::Functions::LocalizableFunctionInterface or a ConstDiscreteFunction.
   */
  template <class SourceType, class RangeSpaceType, class VectorType>
  void apply(const SourceType& source, LocalDiscreteFunction<RangeSpaceType, VectorType>& local_range) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply(source, local_range));
  }
}; // class LocalOperatorInterface


template <class Traits>
class LocalCouplingOperatorInterface : public XT::CRTPInterface<LocalCouplingOperatorInterface<Traits>, Traits>,
                                       public XT::Common::ParametricInterface,
                                       internal::IsLocalCouplingOperator
{
public:
  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_neighbor,
             const XT::Common::Parameter& mu = {}) const
  {
    this->as_imp().apply(source, intersection, local_range_entity, local_range_neighbor, mu);
  }
}; // class LocalCouplingOperatorInterface


template <class Traits>
class LocalBoundaryOperatorInterface : public XT::CRTPInterface<LocalBoundaryOperatorInterface<Traits>, Traits>,
                                       internal::IsLocalBoundaryOperator
{
public:
  template <class SourceType, class IntersectionType, class SpaceType, class VectorType>
  void apply(const SourceType& source,
             const IntersectionType& intersection,
             LocalDiscreteFunction<SpaceType, VectorType>& local_range_entity) const
  {
    this->as_imp().apply(source, intersection, local_range_entity);
  }
}; // class LocalBoundaryOperatorInterface
#endif // 0


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


/// \todo move to type_traits.hh
template <class T>
struct is_local_operator : public std::is_base_of<internal::IsLocalOperator, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_coupling_operator : public std::is_base_of<internal::IsLocalCouplingOperator, T>
{
};

/// \todo move to type_traits.hh
template <class T>
struct is_local_boundary_operator : public std::is_base_of<internal::IsLocalBoundaryOperator, T>
{
};
#endif // 0


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_OPERATORS_INTERFACES_HH
