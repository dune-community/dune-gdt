// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PRODUCTS_PROLONG_HH
#define DUNE_GDT_TEST_OPERATORS_PRODUCTS_PROLONG_HH

#include "operators_products.hh"

// +---------------------------------------+
// |  * to test the prolongation operators |
// +---------------------------------------+

template< class CoarseSpaceType, class FineSpaceType, class ProlongationOperatorType >
struct ProlongationOperatorBase
{
  typedef typename FineSpaceType::GridViewType      GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename FineSpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                       dimDomain = FineSpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename FineSpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                       dimRange = FineSpaceType::dimRange;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 2u);
    auto& grid = grid_provider.grid();
    grid.globalRefine(1);
    const auto coarse_grid_part_view = Dune::GDT::SpaceTools::GridPartView< CoarseSpaceType >::create_level(grid, 0);
    assert(grid.maxLevel() > 0);
    const auto fine_grid_part_view = Dune::GDT::SpaceTools::GridPartView< FineSpaceType >::create_level(grid,
                                                                                                        grid.maxLevel());
    assert(fine_grid_part_view->indexSet().size(0) > coarse_grid_part_view->indexSet().size(0));
    // first, project an anlytical function onto the coarse grid
    const FunctionType function("x", "x[0]", 1, "function");
    const CoarseSpaceType coarse_space(coarse_grid_part_view);
    VectorType coarse_vector(coarse_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< CoarseSpaceType, VectorType > CoarseDiscreteFunctionType;
    CoarseDiscreteFunctionType coarse_discrete_function(coarse_space, coarse_vector, "coarse discrete function");
    const Dune::GDT::Operators::Projection< GridViewType > coarse_projection_operator(*(coarse_space.grid_view()));
    coarse_projection_operator.apply(function, coarse_discrete_function);
    // since the projection operator was tested above we are confident this worked
    // but we check anyway (the L2 product operator was also tested above)
    const Dune::GDT::Products::L2< GridViewType > coarse_l2_product_operator(*(coarse_space.grid_view()));
    const Dune::Stuff::Functions::Difference< FunctionType, CoarseDiscreteFunctionType >
        coarse_difference(function, coarse_discrete_function);
    const auto coarse_l2_error = std::sqrt(coarse_l2_product_operator.apply2(coarse_difference, coarse_difference));
    if (coarse_l2_error > RangeFieldType(1e-15))
      DUNE_THROW(Dune::Stuff::Exceptions::internal_error,
                            "This should not happen, those operators were tested above!\n"
                            << coarse_l2_error << " vs. " << RangeFieldType(1e-15));
    // now we prolong the discrete function from the coarse to the fine grid part
    const FineSpaceType fine_space(fine_grid_part_view);
    VectorType fine_vector(fine_space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< FineSpaceType, VectorType > FineDiscreteFunctionType;
    FineDiscreteFunctionType fine_discrete_function(fine_space, fine_vector, "fine discrete function");
    const ProlongationOperatorType prolongation_operator(*(fine_space.grid_view()));
    prolongation_operator.apply(coarse_discrete_function, fine_discrete_function);
    // and measure the error
    const Dune::GDT::Products::L2< GridViewType > fine_l2_product_operator(*(fine_space.grid_view()));
    const Dune::Stuff::Functions::Difference< FunctionType, FineDiscreteFunctionType >
        fine_difference(function, fine_discrete_function);
    const auto fine_l2_error = std::sqrt(fine_l2_product_operator.apply2(fine_difference, fine_difference));
    EXPECT_LE(fine_l2_error, RangeFieldType(1e-15));
  }
}; // ProlongationOperatorBase

template< class P >
struct L2ProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::L2Prolongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};

template< class P >
struct LagrangeProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::LagrangeProlongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};

template< class P >
struct ProlongationOperator
  : public ProlongationOperatorBase< typename P::first_type,
                                     typename P::second_type,
                                     Dune::GDT::Operators::Prolongation< typename P::second_type::GridViewType > >
  , public ::testing::Test
{};


// +---------------------------------------------------------+
// |  * arguments for the prolongation operator test structs |
// +---------------------------------------------------------+



#define LAGRANGE_PROLONGATION_OPERATOR_SPACE_TYPES \
  /* all combinations which have Spaces::ContinuousLagrange::FemBased as FineSpaceType */ \
    std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLevelGridPartType, 1, double, 1 > > \
  , std::pair< Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLevelGridPartType, 1, double, 1 >, \
               Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLevelGridPartType, 1, double, 1 > >

#endif // DUNE_GDT_TEST_OPERATORS_PRODUCTS_PROLONG_HH
