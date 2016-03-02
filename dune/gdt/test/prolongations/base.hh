// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_BASE_HH
#define DUNE_GDT_TEST_PROLONGATIONS_BASE_HH

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/projections.hh>
#include <dune/gdt/products/l2.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \note Assumes that project and Products::L2 does the right thing.
 */
template< class CoarseSpaceType,
          class FineSpaceType,
          template< class, class, class > class ProlongationOperatorType >
struct LocalizableProlongationOperatorBase
  : public ::testing::Test
{
  typedef typename FineSpaceType::GridViewType      GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const size_t                             dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename FineSpaceType::RangeFieldType  RangeFieldType;
  static const size_t                             dimRange = FineSpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef typename Dune::Stuff::LA::Container< RangeFieldType >::VectorType VectorType;

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15) const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 2u);
    grid_provider.grid().globalRefine(1);
    auto coarse_grid_view = grid_provider.template level< CoarseSpaceType::part_view_type >(0);
    auto fine_grid_view = grid_provider.template level< FineSpaceType::part_view_type >(grid_provider.grid().maxLevel());
    // first, project an anlytical function onto the coarse grid
    const FunctionType function("x", "x[0]", 1, "function");
    const CoarseSpaceType coarse_space(coarse_grid_view);
    typedef Dune::GDT::DiscreteFunction< CoarseSpaceType, VectorType > CoarseDiscreteFunctionType;
    CoarseDiscreteFunctionType coarse_discrete_function(coarse_space, "coarse discrete function");
    project(function, coarse_discrete_function);
    // since the projection operator is tested elsewhere we are confident this worked, but we check anyway
    const Dune::GDT::Products::L2< GridViewType > coarse_l2_product_operator(coarse_space.grid_view());
    const auto coarse_l2_error = coarse_l2_product_operator.induced_norm(function - coarse_discrete_function);
    if (coarse_l2_error > tolerance)
      DUNE_THROW(Dune::Stuff::Exceptions::internal_error,
                            "This should not happen, the L2 product and projection operators are tested elsewhere!\n"
                            << coarse_l2_error << " vs. " << RangeFieldType(1e-15));
    // now we prolong the discrete function from the coarse to the fine grid part
    const FineSpaceType fine_space(fine_grid_view);
    typedef Dune::GDT::DiscreteFunction< FineSpaceType, VectorType > FineDiscreteFunctionType;
    FineDiscreteFunctionType fine_discrete_function(fine_space, "fine discrete function");
    ProlongationOperatorType< GridViewType, CoarseDiscreteFunctionType, FineDiscreteFunctionType >
        prolongation_operator(fine_space.grid_view(), coarse_discrete_function, fine_discrete_function);
    prolongation_operator.apply();
    // and measure the error
    const Dune::GDT::Products::L2< GridViewType > fine_l2_product_operator(fine_space.grid_view());
    const auto fine_l2_error = fine_l2_product_operator.induced_norm(function - fine_discrete_function);
    EXPECT_LE(fine_l2_error, tolerance);
  }
}; // LocalizableProlongationOperatorBase


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_BASE_HH
