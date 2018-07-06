// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_PROJECTIONS_BASE_HH
#define DUNE_GDT_TEST_PROJECTIONS_BASE_HH

#include <dune/common/unused.hh>

#include <dune/xt/common/test/float_cmp.hh>

#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/test/operators/base.hh>

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


/**
 * \note This test assumes that Products::L2 does the right thing!
 */
template <class SpaceType>
struct ProjectionOperatorBase : public OperatorBase<SpaceType>
{
  typedef OperatorBase<SpaceType> BaseType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::GridLayerType;

  void measure_error(const RangeFieldType& expected_error) const
  {
    const auto l2_error =
        make_l2_operator(this->space_.grid_layer(), 2)->induced_norm(this->function_ - this->discrete_function_);
    DXTC_EXPECT_FLOAT_LE(l2_error, expected_error);
  }

  ProjectionOperatorBase()
    : BaseType()
    , relax_factor(this->space_.grid_layer().grid().comm().size() > 1 ? 1.15 : 1)
    , default_tolerance(1.35e-10 * relax_factor)
    , alugrid_tolerance(3.8e-11)
  {
  }

  const double relax_factor;
  const double default_tolerance;
  const double alugrid_tolerance;
}; // struct ProjectionOperatorBase

} // namespace internal


template <class SpaceType, class LocalizableProjectionOperatorType>
struct LocalizableProjectionOperatorBase : public internal::ProjectionOperatorBase<SpaceType>
{
  typedef internal::ProjectionOperatorBase<SpaceType> BaseType;
  using typename BaseType::RangeFieldType;

  void constructible_by_ctor()
  {
    DUNE_UNUSED LocalizableProjectionOperatorType projection_operator(
        this->space_.grid_layer(), this->function_, this->discrete_function_);
  }

  void produces_correct_results(const RangeFieldType& tolerance = BaseType::default_tolerance)
  {
    this->discrete_function_.vector() *= 0.0;
    LocalizableProjectionOperatorType projection_operator(
        this->space_.grid_layer(), this->function_, this->discrete_function_);
    projection_operator.apply();
    this->measure_error(tolerance);
  }
}; // class LocalizableProjectionOperatorBase


template <class SpaceType, class ProjectionOperatorType>
struct ProjectionOperatorBase : public internal::ProjectionOperatorBase<SpaceType>
{
  typedef internal::ProjectionOperatorBase<SpaceType> BaseType;
  using typename BaseType::RangeFieldType;

  void constructible_by_ctor()
  {
    DUNE_UNUSED ProjectionOperatorType projection_operator(this->space_.grid_layer());
  }

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    this->discrete_function_.vector() *= 0.0;
    ProjectionOperatorType projection_operator(this->space_.grid_layer());
    projection_operator.apply(this->function_, this->discrete_function_);
    this->measure_error(tolerance);
  }
}; // class ProjectionOperatorBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_BASE_HH
