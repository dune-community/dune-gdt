// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_PROJECTIONS_BASE_HH
#define DUNE_GDT_TEST_PROJECTIONS_BASE_HH

#include <dune/common/unused.hh>

#include <dune/gdt/operators/l2.hh>
#include <dune/gdt/spaces/tools.hh>

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
  using typename BaseType::GridViewType;

  void measure_error(const RangeFieldType& tolerance) const
  {
    const auto l2_error =
        make_l2_operator(this->space_.grid_view(), 2)->induced_norm(this->function_ - this->discrete_function_);
    EXPECT_LE(l2_error, tolerance) << "l2_error:  " << std::scientific << l2_error << "\n"
                                   << "tolerance: " << std::scientific << tolerance;
  }

  static const constexpr double default_tolerance = 1e-15;
  static const constexpr double alugrid_tolerance = 3.8e-11;
}; // struct ProjectionOperatorBase

template <class T>
constexpr double ProjectionOperatorBase<T>::default_tolerance;
template <class T>
constexpr double ProjectionOperatorBase<T>::alugrid_tolerance;

} // namespace internal


template <class SpaceType, class LocalizableProjectionOperatorType>
struct LocalizableProjectionOperatorBase : public internal::ProjectionOperatorBase<SpaceType>
{
  typedef internal::ProjectionOperatorBase<SpaceType> BaseType;
  using typename BaseType::RangeFieldType;

  void constructible_by_ctor()
  {
    LocalizableProjectionOperatorType DUNE_UNUSED(projection_operator)(
        this->space_.grid_view(), this->function_, this->discrete_function_);
  }

  void produces_correct_results(const RangeFieldType& tolerance = double(BaseType::default_tolerance))
  {
    this->discrete_function_.vector() *= 0.0;
    LocalizableProjectionOperatorType projection_operator(
        this->space_.grid_view(), this->function_, this->discrete_function_);
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
    ProjectionOperatorType DUNE_UNUSED(projection_operator)(this->space_.grid_view());
  }

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    this->discrete_function_.vector() *= 0.0;
    ProjectionOperatorType projection_operator(this->space_.grid_view());
    projection_operator.apply(this->function_, this->discrete_function_);
    this->measure_error(tolerance);
  }
}; // class ProjectionOperatorBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROJECTIONS_BASE_HH
