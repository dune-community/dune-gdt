// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PROJECTION_HH
#define DUNE_GDT_TEST_OPERATORS_PROJECTION_HH

#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/test/gtest/gtest.h>

#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/spaces/tools.hh>

using namespace Dune;
using namespace Dune::GDT;


/**
 * \note This test assumes that Products::L2 does the right thing!
 */
template <class SpaceType, class ProjectionOperatorType>
struct ProjectionOperatorBase : ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Functions::Expression<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange>
      FunctionType;
  typedef typename Stuff::LA::Container<RangeFieldType, Stuff::LA::default_backend>::VectorType VectorType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 3u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const ProjectionOperatorType projection_operator(space.grid_view());
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference<FunctionType, DiscreteFunctionType> difference(function,
                                                                                            discrete_function);
    const Dune::GDT::Products::L2<GridViewType> l2_product_operator(space.grid_view());
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    EXPECT_LE(l2_error, RangeFieldType(1e-11)); // 3d needs this
  }
}; // ProjectionOperatorType


template <class SpaceType>
struct ProjectionOperator
    : public ProjectionOperatorBase<SpaceType, Operators::Projection<typename SpaceType::GridViewType>>
{
  typedef ProjectionOperatorBase<SpaceType, Operators::Projection<typename SpaceType::GridViewType>> BaseType;
  typedef typename BaseType::GridProviderType GridProviderType;
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::FunctionType FunctionType;
  typedef typename BaseType::RangeFieldType RangeFieldType;
  typedef typename BaseType::VectorType VectorType;

  void apply_projection_works() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 3u);
    auto& grid                = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    Operators::apply_projection(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference<FunctionType, DiscreteFunctionType> difference(function,
                                                                                            discrete_function);
    const Dune::GDT::Products::L2<GridViewType> l2_product_operator(space.grid_view());
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    EXPECT_LE(l2_error, RangeFieldType(1e-15));
  }
};


#endif // DUNE_GDT_TEST_OPERATORS_PROJECTION_HH
