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

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/xt/common/crtp.hh>

#include <dune/xt/functions/type_traits.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {


template <class TestBase, class Field = typename TestBase::RangeFieldType>
class LocalVolumeFunctionalInterface
{
  static_assert(XT::Functions::is_localfunction_set<TestBase>::value, "");

public:
  typedef TestBase TestBaseType;
  typedef Field FieldType;

  virtual ~LocalVolumeFunctionalInterface() = default;

  virtual void apply(const TestBaseType& test_basis, DynamicVector<FieldType>& ret) const = 0;

  DynamicVector<FieldType> apply(const TestBaseType& test_basis) const
  {
    DynamicVector<FieldType> ret(test_basis.size(), 0.);
    apply(test_basis, ret);
    return ret;
  }
}; // class LocalFunctionalInterface


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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FUNCTIONALS_INTERFACES_HH
