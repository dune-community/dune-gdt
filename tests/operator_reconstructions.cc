// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>
#include <utility>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
# define ENABLE_ALUGRID 1
# include <dune/grid/alugrid.hh>
#else
# error "This test requires alugrid!"
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/functions/expression.hh>
#include <dune/stuff/functions/combined.hh>

#include <dune/gdt/space/tools.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/playground/operator/darcy.hh>
#include <dune/gdt/product/l2.hh>
#include <dune/gdt/product/h1.hh>

class errors_are_not_as_expected
  : public Dune::Exception
{};

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

template< class SpaceTypes >
struct Darcy_Operator
  : public ::testing::Test
{
#if HAVE_EIGEN
    typedef Dune::Stuff::LA::EigenDenseVector< double > VectorType;
#elif HAVE_DUNE_ISTL
    typedef Dune::Stuff::LA::IstlDenseVector< double>   VectorType;
#else
    typedef Dune::Stuff::LA::CommonDenseVector< double> VectorType;
#endif

  typedef typename SpaceTypes::first_type SourceSpaceType;
  typedef typename SpaceTypes::second_type RangeSpaceType;

  typedef typename RangeSpaceType::GridViewType     GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  static const unsigned int dimDomain = SourceSpaceType::dimDomain;

  void produces_correct_results() const
  {
    using namespace Dune;
    using namespace GDT;

    GridProviderType grid_provider(0.0, 1.0, 4);
    auto grid = grid_provider.grid();
    grid->globalRefine(1);

    typedef Stuff::Function::Expression< typename GridViewType::template Codim< 0 >::Entity
                                       , typename GridViewType::ctype, GridViewType::dimension
                                       , double, 1, 1 > FunctionType;
    const FunctionType source("x", "x[0] * x[1]", 2, "source", {{"x[1]", "x[0]"}});

    const RangeSpaceType range_space(SpaceTools::GridPartView< RangeSpaceType >::create_leaf(*grid));
    VectorType range_vector(range_space.mapper().size());
    DiscreteFunction< RangeSpaceType, VectorType > range(range_space, range_vector);

    const FunctionType function("x", "-1.0", 0);
    const GDT::Operator::DarcyReconstruction< GridViewType, FunctionType > darcy_operator(*(range_space.grid_view()),
                                                                                          function);
    darcy_operator.apply(source, range);

    const Stuff::Function::Expression< typename GridViewType::template Codim< 0 >::Entity
                                     , typename GridViewType::ctype, GridViewType::dimension
                                     , double, 2, 1 >
      desired_output("x", std::vector< std::string >({"x[1]", "x[0]"}), 1,
                     "desired output",
                     {{"0.0", "1.0"}, {"1.0", "0.0"}});
    const auto difference = desired_output - range;

    const Product::L2Generic< GridViewType > l2_product(*(range_space.grid_view()));
    const double l2_error = std::sqrt(l2_product.apply2(difference, difference));
    if (l2_error > 1e-15)
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected, l2_error);

    const Product::H1SemiGeneric< GridViewType > h1_semi_product(*(range_space.grid_view()));
    const double h1_semi_error = std::sqrt(h1_semi_product.apply2(difference, difference));
    if (h1_semi_error > 1e-14)
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected, h1_semi_error);
  } // ... produces_correct_results()
}; // struct Darcy_Operator

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

typedef Dune::ALUConformGrid< 2, 2 > AluConform2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, false >::Type AluConform2dLeafGridPartType;

typedef testing::Types< std::pair< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dLeafGridPartType, 1, double, 1 >,
                                   Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dLeafGridPartType, 1, double, 2 > >
                      > SpaceTypes;

// +--------------------------------------------------------------------------------------+
// | 3rd we combine all test structs with their appropriate arguments to create the tests |
// | (comment out the following lines if you do not want a test to be run)                |
// +--------------------------------------------------------------------------------------+

TYPED_TEST_CASE(Darcy_Operator, SpaceTypes);
TYPED_TEST(Darcy_Operator, produces_correct_results) {
  this->produces_correct_results();
}

// +--------------------------------------------------------------------------------------+
// | 4th we run all the tests                                                             |
// | (run the resulting executable with '--gtest_catch_exceptions=0' to see an exception) |
// +--------------------------------------------------------------------------------------+

int main(int argc, char** argv)
{
  try {
    test_init(argc, argv);
    return RUN_ALL_TESTS();
  } catch (Dune::Exception& e) {
    std::cerr << "\nDune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << "\n" << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // ... main(...)
