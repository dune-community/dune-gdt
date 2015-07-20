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


/**
 * \note This test assumes that Products::L2 does the right thing!
 */
template< class SpaceType, class ProjectionOperatorType >
class ProjectionOperatorBase
  : public ::testing::Test
{
protected:
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const size_t                         dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int polOrder = SpaceType::polOrder;
  static const size_t                         dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;
  typedef typename Dune::Stuff::LA::Container< RangeFieldType, Dune::Stuff::LA::default_backend >::VectorType VectorType;

public:
  ProjectionOperatorBase()
    : grid_provider_(0.0, 1.0, 3u)
    , space_(Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid_provider_.grid()))
    , function_("x", {"x[0]", "0", "0"}, 1)
    , vector_(space_.mapper().size())
    , discrete_function_(space_, vector_)
  {}

  void produces_correct_results(const RangeFieldType& tolerance = 1e-15)
  {
    vector_ *= 0.0;
    const ProjectionOperatorType projection_operator(space_.grid_view());
    projection_operator.apply(function_, discrete_function_);
    measure_error(tolerance);
  }

  void free_project_function_works(const RangeFieldType& tolerance = 1e-15)
  {
    vector_ *= 0.0;
    Dune::GDT::project(function_, discrete_function_);
    measure_error(tolerance);
  }

protected:
  void measure_error(const RangeFieldType& tolerance) const
  {
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(space_.grid_view());
    const auto l2_error = l2_product_operator.induced_norm(function_ - discrete_function_);
    EXPECT_LE(l2_error, tolerance);
  }

  GridProviderType grid_provider_;
  const SpaceType space_;
  const FunctionType function_;
  VectorType vector_;
  Dune::GDT::DiscreteFunction< SpaceType, VectorType > discrete_function_;
}; // class ProjectionOperatorBase


template< class SpaceType >
struct ProjectionOperator
  : ProjectionOperatorBase< SpaceType, Dune::GDT::Operators::Projection< typename SpaceType::GridViewType > >
{};


#endif // DUNE_GDT_TEST_OPERATORS_PROJECTION_HH
