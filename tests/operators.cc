// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>

#include <dune/common/exceptions.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>
#include <dune/stuff/la/container/eigen.hh>

#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operator/projections.hh>
#include <dune/gdt/operator/products.hh>

class errors_are_not_as_expected
  : public Dune::Exception
{};

typedef Dune::Stuff::LA::EigenDenseVector< double > VectorType;

typedef Dune::SGrid< 1, 1 >                           S1dGridType;
typedef Dune::grid::Part::Leaf::Const< S1dGridType >  S1dGridPartType;
typedef Dune::SGrid< 2, 2 >                           S2dGridType;
typedef Dune::grid::Part::Leaf::Const< S2dGridType >  S2dGridPartType;
typedef Dune::SGrid< 3, 3 >                           S3dGridType;
typedef Dune::grid::Part::Leaf::Const< S3dGridType >  S3dGridPartType;

typedef Dune::YaspGrid< 1 >                             Yasp1dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp1dGridType > Yasp1dGridPartType;
typedef Dune::YaspGrid< 2 >                             Yasp2dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp2dGridType > Yasp2dGridPartType;
typedef Dune::YaspGrid< 3 >                             Yasp3dGridType;
typedef Dune::grid::Part::Leaf::Const< Yasp3dGridType > Yasp3dGridPartType;

#if HAVE_ALUGRID
typedef Dune::ALUConformGrid< 2, 2 >                          AluConform2dGridType;
typedef Dune::grid::Part::Leaf::Const< AluConform2dGridType > AluConform2dGridPartType;
typedef Dune::ALUSimplexGrid< 2, 2 >                          AluSimplex2dGridType;
typedef Dune::grid::Part::Leaf::Const< AluSimplex2dGridType > AluSimplex2dGridPartType;
typedef Dune::ALUSimplexGrid< 3, 3 >                          AluSimplex3dGridType;
typedef Dune::grid::Part::Leaf::Const< AluSimplex3dGridType > AluSimplex3dGridPartType;
typedef Dune::ALUCubeGrid< 3, 3 >                             AluCube3dGridType;
typedef Dune::grid::Part::Leaf::Const< AluCube3dGridType >    AluCube3dGridPartType;
#endif


typedef testing::Types< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex3dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluCube3dGridPartType, 1, double, 1 >
#endif
                      > L2ProductOperatorSpaceTypes;

template< class SpaceType >
struct L2ProductOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const Dune::GDT::ProductOperator::L2< GridPartType > l2_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    function_1.visualize(grid_part->gridView(), "function1");
    auto l2_product = l2_product_operator.apply2(function_1, function_1);
    RangeFieldType error = l2_product - RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0)
                 << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    l2_product = l2_product_operator.apply2(function_2, function_2);
    error = l2_product - RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/3.0)
                 << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    l2_product = l2_product_operator.apply2(function_3, function_3);
    error = l2_product - RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_product << " vs. " << RangeFieldType(1.0/5.0)
                 << " (difference: " << std::scientific << error << ")");
  }
}; // L2ProductOperator

TYPED_TEST_CASE(L2ProductOperator, L2ProductOperatorSpaceTypes);
TYPED_TEST(L2ProductOperator, produces_correct_results) {
  this->produces_correct_results();
}


template< class SpaceType >
struct H1SemiProductOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  typedef Dune::Stuff::Function::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const Dune::GDT::ProductOperator::H1Semi< GridPartType > h1semi_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    auto h1semi_product = h1semi_product_operator.apply2(function_1, function_1);
    RangeFieldType error = h1semi_product - dimDomain * RangeFieldType(1.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain *RangeFieldType(1.0)
                 << " (difference: " << std::scientific << error << ")");
    // test 2 (linear)
    const FunctionType function_2("x", "fake_value", 2, "affine gradient",
                                  {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    h1semi_product = h1semi_product_operator.apply2(function_2, function_2);
    error = h1semi_product - dimDomain * RangeFieldType(1.0/3.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain * RangeFieldType(1.0/3.0)
                 << " (difference: " << std::scientific << error << ")");
    // test 3 (quadratic)
    const FunctionType function_3("x", "fake_value", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    h1semi_product = h1semi_product_operator.apply2(function_3, function_3);
    error = h1semi_product - dimDomain * RangeFieldType(1.0/5.0);
    if (error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << h1semi_product << " vs. " << dimDomain * RangeFieldType(1.0/5.0)
                 << " (difference: " << std::scientific << error << ")");
  }
}; // H1SemiProductOperator

TYPED_TEST_CASE(H1SemiProductOperator, L2ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiProductOperator, produces_correct_results) {
  this->produces_correct_results();
}


typedef testing::Types< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex3dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluCube3dGridPartType, 1, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 1, double, 1 >
#endif
                      > LagrangeProjectionOperatorSpaceTypes;

template< class SpaceType >
struct LagrangeProjectionOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Function::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::ProjectionOperator::Lagrange< GridPartType > lagrange_projection_operator(*grid_part);
    lagrange_projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::ProductOperator::L2< GridPartType > l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // L2ProjectionOperator

TYPED_TEST_CASE(LagrangeProjectionOperator, LagrangeProjectionOperatorSpaceTypes);
TYPED_TEST(LagrangeProjectionOperator, produces_correct_results) {
  this->produces_correct_results();
}


typedef testing::Types<
#if HAVE_ALUGRID
                        Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 1, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 2, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 2, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 2, double, 1 >
#endif
                      > L2ProjectionOperatorSpaceTypes;

template< class SpaceType >
struct L2ProjectionOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Function::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::ProjectionOperator::L2< GridPartType > l2_projection_operator(*grid_part);
    l2_projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::ProductOperator::L2< GridPartType > l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // L2ProjectionOperator

TYPED_TEST_CASE(L2ProjectionOperator, L2ProjectionOperatorSpaceTypes);
TYPED_TEST(L2ProjectionOperator, produces_correct_results) {
  this->produces_correct_results();
}


typedef testing::Types< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluSimplex3dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluCube3dGridPartType, 1, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluSimplex3dGridPartType, 1, double, 1 >
#endif
                      > DirichletProjectionOperatorSpaceTypes;

template< class SpaceType >
struct DirichletProjectionOperator
  : public ::testing::Test
{
  typedef typename SpaceType::GridPartType          GridPartType;
  typedef typename GridPartType::GridType           GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int polOrder = SpaceType::polOrder;
  typedef Dune::Stuff::Function::Expression
      < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange > FunctionType;

  void produces_correct_results() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 1u); // this has to be 1, otherwise the projection does not equal
    const auto grid = grid_provider.grid();             // x[0] any more!
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    DomainType dirichlet_normal(0);
    dirichlet_normal[0] = DomainFieldType(1);
    const Dune::Stuff::GridboundaryNormalBased< typename GridPartType::IntersectionType >
        boundary_info(false, {dirichlet_normal});
    const SpaceType space(grid_part);
    const FunctionType function("x", "x[0]", 1, "function");
    VectorType vector(space.mapper().size());
    typedef Dune::GDT::DiscreteFunction< SpaceType, VectorType > DiscreteFunctionType;
    DiscreteFunctionType discrete_function(space, vector, "discrete function");
    // project
    const Dune::GDT::ProjectionOperator::Dirichlet< GridPartType > l2_projection_operator(*grid_part, boundary_info);
    l2_projection_operator.apply(function, discrete_function);
    // measure error
    const Dune::Stuff::Function::Difference< FunctionType, DiscreteFunctionType > difference(function,
                                                                                             discrete_function);
    const Dune::GDT::ProductOperator::L2< GridPartType > l2_product_operator(*grid_part);
    const auto l2_error = std::sqrt(l2_product_operator.apply2(difference, difference));
    if (l2_error > RangeFieldType(1e-15))
      DUNE_THROW(errors_are_not_as_expected,
                 "They really ain't!\n" << l2_error << " vs. " << RangeFieldType(1e-15));
  }
}; // DirichletProjectionOperator

TYPED_TEST_CASE(DirichletProjectionOperator, DirichletProjectionOperatorSpaceTypes);
TYPED_TEST(DirichletProjectionOperator, produces_correct_results) {
  this->produces_correct_results();
}


int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
}
