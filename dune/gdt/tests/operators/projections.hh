// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PROJECTIONS_HH
#define DUNE_GDT_TEST_OPERATORS_PROJECTIONS_HH

#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/tools.hh>

#include "../operators.hh"

namespace Dune {
namespace GDT {
namespace Tests {


/**
 * \note This test assumes that Products::L2 does the right thing!
 */
template< class SpaceType, class ProjectionOperatorType >
struct ProjectionOperatorsBase
  : public OperatorBase< SpaceType >
{
  typedef OperatorBase< SpaceType > BaseType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::GridViewType;

  void measure_error(const RangeFieldType& tolerance) const
  {
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(this->space_.grid_view());
    const auto l2_error = l2_product_operator.induced_norm(this->function_ - this->discrete_function_);
    EXPECT_LE(l2_error, tolerance);
  }
}; // struct ProjectionOperatorsBase


template< class SpaceType, class LocalizableProjectionOperatorType >
class LocalizableProjectionOperatorBase
  : public ProjectionOperatorsBase< SpaceType, LocalizableProjectionOperatorType >
{
  typedef ProjectionOperatorsBase< SpaceType, LocalizableProjectionOperatorType > BaseType;
public:
  using typename BaseType::RangeFieldType;

  void constructible_by_ctor()
  {
    LocalizableProjectionOperatorType projection_operator(this->space_.grid_view(),
                                                          this->function_,
                                                          this->discrete_function_);
  }

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    this->discrete_function_.vector() *= 0.0;
    LocalizableProjectionOperatorType projection_operator(this->space_.grid_view(),
                                                          this->function_,
                                                          this->discrete_function_);
    projection_operator.apply();
    this->measure_error(tolerance);
  }
}; // class LocalizableProjectionOperatorBase


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_PROJECTIONS_HH
