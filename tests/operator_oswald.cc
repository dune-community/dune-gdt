// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

// This one has to come first (includes the config.h)!
#include <dune/stuff/test/test_common.hh>

#include <memory>

#if HAVE_ALUGRID_SERIAL_H || HAVE_ALUGRID_PARALLEL_H
#define ENABLE_ALUGRID 1
#include <dune/grid/alugrid.hh>
#else
#error "This test requires alugrid!"
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/spaces/discontinuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>

class errors_are_not_as_expected : public Dune::Exception
{
};

typedef Dune::Stuff::LA::CommonDenseVector<double> VectorType;

// +----------------------------------------------------------------------------+
// | 1st we define all the test structs that do something at the end of the day |
// +----------------------------------------------------------------------------+

template <class SpaceType>
struct Oswald_Interpolation_Operator : public ::testing::Test
{
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  //  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  //  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = SpaceType::dimDomain;
  //  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  //  static const unsigned int                   dimRange = SpaceType::dimRange;

  void produces_correct_results() const
  {
    using namespace Dune;
    using namespace GDT;
    // prepare
    const size_t num_partitions = 2;
    GridProviderType grid_provider(0.0, 1.0, num_partitions);
    auto grid = grid_provider.grid();
    grid->globalRefine(1);
    const auto grid_part_view = SpaceTools::GridPartView<SpaceType>::create_leaf(*grid);
    const SpaceType space(grid_part_view);
    typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    VectorType source_vector(space.mapper().size());
    DiscreteFunctionType source(space, source_vector);
    for (auto entity_ptr = grid_part_view->template begin<0>(); entity_ptr != grid_part_view->template end<0>();
         ++entity_ptr) {
      const auto& entity = *entity_ptr;
      const auto center  = entity.geometry().center();
      double value = 0.0;
      if (center[0] > 0.5)
        value                      = 1.0;
      auto local_source            = source.local_discrete_function(entity);
      auto local_source_DoF_vector = local_source.vector();
      for (size_t local_DoF = 0; local_DoF < local_source_DoF_vector.size(); ++local_DoF)
        local_source_DoF_vector.set(local_DoF, value);
    }
    source.visualize("source", false);
    VectorType range_vector(space.mapper().size());
    DiscreteFunctionType range(space, range_vector);
    Operators::OswaldInterpolation<typename SpaceType::GridViewType> oswald_operator(*(space.grid_view()));
    oswald_operator.apply(source, range);
    range.visualize("range", false);
  } // ... produces_correct_results()
}; // struct Oswald_Interpolation_Operator

// +----------------------------------------------------------------------------+
// | 2nd we define all arguments the above test structs are to be compiled with |
// +----------------------------------------------------------------------------+

typedef Dune::ALUConformGrid<2, 2> AluConform2dGridType;
typedef
    typename Dune::GDT::SpaceTools::LeafGridPartView<AluConform2dGridType, false>::Type AluConform2dLeafGridPartType;
// typedef Dune::ALUSimplexGrid< 2, 2 > AluSimplex2dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex2dGridType, false >::Type
// AluSimplex2dLeafGridPartType;
// typedef Dune::ALUSimplexGrid< 3, 3 > AluSimplex3dGridType;
// typedef typename Dune::GDT::SpaceTools::LeafGridPartView< AluSimplex3dGridType, false >::Type
// AluSimplex3dLeafGridPartType;

typedef testing::Types<Dune::GDT::Spaces::DiscontinuousLagrange::FemLocalfunctionsBased<AluConform2dLeafGridPartType, 1,
                                                                                        double, 1>> SpaceTypes;

// +--------------------------------------------------------------------------------------+
// | 3rd we combine all test structs with their appropriate arguments to create the tests |
// | (comment out the following lines if you do not want a test to be run)                |
// +--------------------------------------------------------------------------------------+


TYPED_TEST_CASE(Oswald_Interpolation_Operator, SpaceTypes);
TYPED_TEST(Oswald_Interpolation_Operator, produces_correct_results)
{
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
}
