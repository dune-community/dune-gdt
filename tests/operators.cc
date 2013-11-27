// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// disable all grids before we include config.h
#ifdef GRIDDIM
# undef GRIDDIM
#endif
#ifdef YASPGRID
# undef YASPGRID
#endif
#ifdef ALBERTAGRID
# undef ALBERTAGRID
#endif
#ifdef UGGRID
# undef UGGRID
#endif
#ifdef ALUGRID_CONFORM
# undef ALUGRID_CONFORM
#endif
#ifdef ALUGRID_CUBE
# undef ALUGRID_CUBE
#endif
#ifdef ALUGRID_SIMPLEX
# undef ALUGRID_SIMPLEX
#endif
#ifdef ONEDGRID
# undef ONEDGRID
#endif
#ifdef SGRID
# undef SGRID
#endif

#include <dune/stuff/test/test_common.hh>

# include <memory>

# include <dune/common/exceptions.hh>
# include <dune/common/float_cmp.hh>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
// enable alugrid
# define GRIDDIM 2
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#endif
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <dune/grid/part/leaf.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
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

typedef Dune::ALUConformGrid< 2, 2 >                                  AluConform2dGridType;
typedef Dune::grid::Part::Leaf::Const< AluConform2dGridType >         AluConform2dGridPartType;


typedef testing::Types< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S1dGridPartType, 2, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S2dGridPartType, 2, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< S3dGridPartType, 2, double, 1 >

                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp1dGridPartType, 2, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp2dGridPartType, 2, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< Yasp3dGridPartType, 2, double, 1 >
#if HAVE_ALUGRID
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dGridPartType, 2, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::ContinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 2, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 1, double, 1 >
                      , Dune::GDT::DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< AluConform2dGridPartType, 2, double, 1 >
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

  void check() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const Dune::GDT::ProductOperator::L2< GridPartType > l2_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "1.0", 0);
    auto l2_product = l2_product_operator.apply2(function_1, function_1);
    if (Dune::FloatCmp::ne(l2_product, RangeFieldType(1)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
    // test 2 (linear)
    const FunctionType function_2("x", "x[0] - 1.0", 1);
    l2_product = l2_product_operator.apply2(function_2, function_2);
    if (Dune::FloatCmp::ne(l2_product, RangeFieldType(1.0/3.0)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
    // test 3 (quadratic)
    const FunctionType function_3("x", "x[0]*x[0]", 2);
    l2_product = l2_product_operator.apply2(function_3, function_3);
    if (Dune::FloatCmp::ne(l2_product, RangeFieldType(1.0/5.0)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
  }
};

TYPED_TEST_CASE(L2ProductOperator, L2ProductOperatorSpaceTypes);
TYPED_TEST(L2ProductOperator, produces_correct_results) {
  this->check();
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

  void check() const
  {
    // prepare
    const GridProviderType grid_provider(0.0, 1.0, 4u);
    const auto grid = grid_provider.grid();
    const auto grid_part = std::make_shared< const GridPartType >(*grid);
    const Dune::GDT::ProductOperator::H1Semi< GridPartType > h1semi_product_operator(*grid_part);
    // test 1 (constant)
    const FunctionType function_1("x", "fake_value", 1, "constant gradient", {{"1.0", "1.0", "1.0"}});
    auto h1semi_product = h1semi_product_operator.apply2(function_1, function_1);
    if (Dune::FloatCmp::ne(h1semi_product, dimDomain * RangeFieldType(1.0)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
    // test 2 (linear)
    const FunctionType function_2("x", "fake_value", 2, "affine gradient",
                                  {{"x[0] - 1.0", "x[0] - 1.0", "x[0] - 1.0"}});
    h1semi_product = h1semi_product_operator.apply2(function_2, function_2);
    if (Dune::FloatCmp::ne(h1semi_product, dimDomain * RangeFieldType(1.0/3.0)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
    // test 3 (quadratic)
    const FunctionType function_3("x", "fake_value", 3, ", quadratic gradient",
                                  {{"x[0]*x[0]", "x[0]*x[0]", "x[0]*x[0]"}});
    h1semi_product = h1semi_product_operator.apply2(function_3, function_3);
    if (Dune::FloatCmp::ne(h1semi_product, dimDomain * RangeFieldType(1.0/5.0)))
      DUNE_THROW(errors_are_not_as_expected, "They really ain't!");
  }
};

TYPED_TEST_CASE(H1SemiProductOperator, L2ProductOperatorSpaceTypes);
TYPED_TEST(H1SemiProductOperator, produces_correct_results) {
  this->check();
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
