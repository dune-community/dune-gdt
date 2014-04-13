// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
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

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/continuouslagrange/fem.hh>
#include <dune/gdt/playground/spaces/finitevolume.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/darcy.hh>
#include <dune/gdt/operators/projections.hh>
#include <dune/gdt/products/l2.hh>
#include <dune/gdt/products/h1.hh>

class errors_are_not_as_expected
  : public Dune::Exception
{};

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

template< class SpaceTypes >
class Darcy_Operator
  : public ::testing::Test
{
  typedef typename SpaceTypes::first_type SourceSpaceType;
  typedef typename SpaceTypes::second_type RangeSpaceType;

  typedef typename RangeSpaceType::GridViewType     GridViewType;
  typedef typename GridViewType::Grid               GridType;
  typedef Dune::Stuff::GridProviderCube< GridType > GridProviderType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = SourceSpaceType::dimDomain;
  typedef double RangeFieldType;

#if HAVE_EIGEN
    typedef Dune::Stuff::LA::EigenDenseVector< RangeFieldType >   VectorType;
#elif HAVE_DUNE_ISTL
    typedef Dune::Stuff::LA::IstlDenseVector< RangeFieldType >    VectorType;
#else
    typedef Dune::Stuff::LA::CommonDenseVector< RangeFieldType >  VectorType;
#endif

public:
  void produces_correct_results() const
  {
    using namespace Dune;
    using namespace GDT;

    GridProviderType grid_provider(0.0, 1.0, 4);
    auto grid = grid_provider.grid();
    grid->globalRefine(1);

    typedef Stuff::Function::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > FunctionType;
    const FunctionType source("x", "x[0] * x[1]", 2, "source", {{"x[1]", "x[0]"}});

    const RangeSpaceType range_space(SpaceTools::GridPartView< RangeSpaceType >::create_leaf(*grid));
    VectorType range_vector(range_space.mapper().size());
    DiscreteFunction< RangeSpaceType, VectorType > range(range_space, range_vector);

    const FunctionType function("x", "-1.0", 0);
    const Operators::Darcy< GridViewType, FunctionType > darcy_operator(*(range_space.grid_view()),
                                                                                          function);
    darcy_operator.apply(source, range);

    const Stuff::Function::Expression< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimDomain >
      desired_output("x", std::vector< std::string >({"x[1]", "x[0]"}), 1,
                     "desired output",
                     {{"0.0", "1.0"}, {"1.0", "0.0"}});

    const Products::L2< GridViewType > l2_product(*(range_space.grid_view()));
    const RangeFieldType l2_error = l2_product.induced_norm(desired_output - range);
    const RangeFieldType l2_error_expected = expected_result_("l2", desired_output, range_space.grid_view());
    if (l2_error > l2_error_expected)
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected, l2_error << " vs. " << l2_error_expected);

    const Products::H1SemiGeneric< GridViewType > h1_semi_product(*(range_space.grid_view()));
    const RangeFieldType h1_error = h1_semi_product.induced_norm(desired_output - range);
    const RangeFieldType h1_error_expected = expected_result_("h1", desired_output, range_space.grid_view());
    if (h1_error > h1_error_expected)
      DUNE_THROW_COLORFULLY(errors_are_not_as_expected, h1_error << " vs. " << h1_error_expected);
  } // ... produces_correct_results()

private:
  template< class FunctionType, class GV >
  RangeFieldType expected_result_(const std::string type,
                                  const FunctionType& desired_output,
                                  const std::shared_ptr< const GV >& grid_view_ptr) const
  {
    typedef typename Dune::GDT::SpaceTools::LeafGridPartView< GridType, RangeSpaceType::needs_grid_view >::Type GPV;
    if (std::is_base_of< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< GPV, 1, RangeFieldType, dimDomain >
                       , RangeSpaceType >::value) {
      if (type == "l2")
        return 2.18e-16;
      else if (type == "h1")
        return 3.12e-15;
      else
        DUNE_THROW_COLORFULLY(Dune::Stuff::Exceptions::internal_error, type);
    } else if (std::is_base_of< Dune::GDT::RaviartThomasSpace::PdelabBased< GPV, 0, RangeFieldType, dimDomain >
                              , RangeSpaceType >::value) {
      typedef Dune::GDT::FiniteVolumeSpace::Default< GV, RangeFieldType, dimDomain > FvSpaceType;
      const FvSpaceType fv_space(grid_view_ptr);
      VectorType fv_desired_output_vector(fv_space.mapper().size());
      Dune::GDT::DiscreteFunction< FvSpaceType, VectorType > fv_desired_output(fv_space, fv_desired_output_vector);
      const Dune::GDT::Operators::L2Projection< GV > l2_projection(*grid_view_ptr);
      l2_projection.apply(desired_output, fv_desired_output);
      const Dune::GDT::Products::L2< GV > l2_product(*grid_view_ptr);
      const Dune::GDT::Products::H1SemiGeneric< GV > h1_semi_product(*grid_view_ptr);
      if (type == "l2")
        return 2.0 * l2_product.induced_norm(desired_output - fv_desired_output);
      else if (type == "h1")
        return h1_semi_product.induced_norm(desired_output - fv_desired_output);
      else
        DUNE_THROW_COLORFULLY(Dune::Stuff::Exceptions::internal_error, type);
    } else
      DUNE_THROW_COLORFULLY(Dune::Stuff::Exceptions::internal_error, type);
  } // ... expected_result_(...)
}; // class Darcy_Operator

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

typedef Dune::ALUConformGrid< 2, 2 > AluConform2dGridType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, true >::Type  AluConform2dLeafGridViewType;
typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluConform2dGridType, false >::Type AluConform2dLeafGridPartType;

typedef testing::Types< std::pair< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dLeafGridPartType, 1, double, 1 >,
                                   Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dLeafGridPartType, 1, double, 2 > >
                      , std::pair< Dune::GDT::ContinuousLagrangeSpace::FemWrapper< AluConform2dLeafGridPartType, 1, double, 1 >,
                                   Dune::GDT::RaviartThomasSpace::PdelabBased< AluConform2dLeafGridViewType, 0, double, 2 > >
                      > SpaceTypes;

// +--------------------------------------------------------------------------------------+
// | 3rd we combine all test structs with their appropriate arguments to create the tests |
// | (comment out the following lines if you do not want a test to run)                   |
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
