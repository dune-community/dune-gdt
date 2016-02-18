// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PROJECTIONS_BASE_HH
#define DUNE_GDT_TEST_OPERATORS_PROJECTIONS_BASE_HH

#include <dune/common/unused.hh>

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/tools.hh>

#include "../base.hh"

namespace Dune {
namespace GDT {
namespace Test {
namespace internal {


/**
 * \note This test assumes that Products::L2 does the right thing!
 */
template <class SpaceType, class ProjectionOperatorType>
struct ProjectionOperatorBase : public OperatorBase<SpaceType>
{
  typedef OperatorBase<SpaceType> BaseType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::GridViewType;

  void measure_error(const RangeFieldType& tolerance) const
  {
    const Dune::GDT::Products::L2<GridViewType> l2_product_operator(this->space_.grid_view());
    const auto l2_error = l2_product_operator.induced_norm(this->function_ - this->discrete_function_);
    EXPECT_LE(l2_error, tolerance);
  }
}; // struct ProjectionOperatorBase


} // namespace internal


template <class SpaceType, class LocalizableProjectionOperatorType>
struct LocalizableProjectionOperatorBase
    : public internal::ProjectionOperatorBase<SpaceType, LocalizableProjectionOperatorType>
{
  typedef internal::ProjectionOperatorBase<SpaceType, LocalizableProjectionOperatorType> BaseType;
  using typename BaseType::RangeFieldType;

  void constructible_by_ctor()
  {
    LocalizableProjectionOperatorType DUNE_UNUSED(projection_operator)(
        this->space_.grid_view(), this->function_, this->discrete_function_);
  }

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    this->discrete_function_.vector() *= 0.0;
    LocalizableProjectionOperatorType projection_operator(
        this->space_.grid_view(), this->function_, this->discrete_function_);
    projection_operator.apply();
    this->measure_error(tolerance);
  }
}; // class LocalizableProjectionOperatorBase


template <class SpaceType, class ProjectionOperatorType>
struct ProjectionOperatorBase : public internal::ProjectionOperatorBase<SpaceType, ProjectionOperatorType>
{
  typedef internal::ProjectionOperatorBase<SpaceType, ProjectionOperatorType> BaseType;
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

#endif // DUNE_GDT_TEST_OPERATORS_PROJECTIONS_BASE_HH
