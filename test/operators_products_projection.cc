// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "operators_products.hh"

// +-------------------------------------+
// |  * to test the projection operators |
// +-------------------------------------+

template< class SpaceType, class ProjectionOperatorType >
struct ProjectionOperatorBase
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 4u);
    auto& grid = grid_provider.grid();
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid);
    const SpaceType space(grid_part_view);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const ProjectionOperatorType projection_operator(*(space.grid_view()));
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(*(space.grid_view()));
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
    Dune::GDT::Operators::apply_projection(function, discrete_function);
  }
}; // ProjectionOperatorType

template< class SpaceType >
struct L2ProjectionOperator
  : public ProjectionOperatorBase< SpaceType, Dune::GDT::Operators::L2Projection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct LagrangeProjectionOperator
  : public ProjectionOperatorBase< SpaceType,
                                   Dune::GDT::Operators::LagrangeProjection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct ProjectionOperator
  : public ProjectionOperatorBase< SpaceType,
                                   Dune::GDT::Operators::Projection< typename SpaceType::GridViewType > >
  , public ::testing::Test
{};

template< class SpaceType >
struct DirichletProjectionOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridViewType          GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::Grid::Providers::Cube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Functions::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    GridProviderType grid_provider(0.0, 1.0, 1u); // this has to be 1, otherwise the projection does not equal
    auto& grid = grid_provider.grid();             // x[0] any more!
    const auto grid_part_view = Dune::GDT::SpaceTools::GridPartView< SpaceType >::create_leaf(grid);
    const SpaceType space(grid_part_view);
    DomainType dirichlet_normal(0);
    dirichlet_normal[0] = DomainFieldType(1);
    const Dune::Stuff::Grid::BoundaryInfos::NormalBased< typename GridViewType::Intersection >
        boundary_info(false, {dirichlet_normal});
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::Operators::DirichletProjection< GridViewType > projection_operator(*(space.grid_view()),
                                                                                       boundary_info);
    projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Functions::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::Products::L2< GridViewType > l2_product_operator(*(space.grid_view()));
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected,
                            "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
    Dune::GDT::Operators::apply_dirichlet_projection(boundary_info, function, discrete_function);
  }
}; // DirichletProjectionOperator


// +-------------------------------------------------------+
// |  * arguments for the projection operator test structs |
// +-------------------------------------------------------+

#define L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID \
    Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluConform2dLeafGridPartType, 2, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 2, double, 1 > \
  , Dune::GDT::Spaces::DiscontinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 2, double, 1 >

typedef testing::Types<
#if HAVE_ALUGRID
                        L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > L2ProjectionOperatorSpaceTypes;

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES \
    Dune::GDT::Spaces::ContinuousLagrange::FemBased< S1dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< S3dLeafGridPartType, 1, double, 1 > \
  \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp1dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< Yasp3dLeafGridPartType, 1, double, 1 >

#define LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID \
    Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluCube3dLeafGridPartType, 1, double, 1 > \
  \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluConform2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex2dLeafGridPartType, 1, double, 1 > \
  , Dune::GDT::Spaces::ContinuousLagrange::FemBased< AluSimplex3dLeafGridPartType, 1, double, 1 >


typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > LagrangeProjectionOperatorSpaceTypes;

typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
                      , L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > ProjectionOperatorSpaceTypes;

typedef testing::Types<
                        LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#if HAVE_ALUGRID
                      , LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#endif
                      > DirichletProjectionOperatorSpaceTypes;

#undef L2_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES
#undef LAGRANGE_PROJECTION_OPERATOR_SPACE_TYPES_ALUGRID


TYPED_TEST_CASE(L2ProjectionOperator, L2ProjectionOperatorSpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

TYPED_TEST_CASE(LagrangeProjectionOperator, LagrangeProjectionOperatorSpaceTypes);
TYPED_TEST(LagrangeProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

TYPED_TEST_CASE(ProjectionOperator, ProjectionOperatorSpaceTypes);
TYPED_TEST(ProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}

TYPED_TEST_CASE(DirichletProjectionOperator, DirichletProjectionOperatorSpaceTypes);
TYPED_TEST(DirichletProjectionOperator, produces_correct_results) {
 this->produces_correct_results();
}


#include <dune/stuff/test/test_main.hh>
