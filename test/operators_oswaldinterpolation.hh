// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_OSWALDINTERPOLATION_HH
#define DUNE_GDT_TEST_OPERATORS_OSWALDINTERPOLATION_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/la/container/common.hh>

#include <dune/gdt/spaces/tools.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>

using namespace Dune;
using namespace GDT;


template <class SpaceType>
struct Oswald_Interpolation_Operator : public ::testing::Test
{
  typedef Dune::Stuff::LA::CommonDenseVector<double> VectorType;
  typedef typename SpaceType::GridViewType GridViewType;
  typedef typename GridViewType::Grid GridType;
  typedef Dune::Stuff::Grid::Providers::Cube<GridType> GridProviderType;
  static const size_t dimDomain = SpaceType::dimDomain;

  /**
   * \note This is not a test for correct results yet!
   */
  void produces_correct_results() const
  {
    // prepare the source function
    const size_t num_partitions = 2;
    GridProviderType grid_provider(0.0, 1.0, num_partitions);
    auto& grid = grid_provider.grid();
    grid.globalRefine(1);
    const auto grid_part_view = SpaceTools::GridPartView<SpaceType>::create_leaf(grid);
    const SpaceType space(grid_part_view);
    typedef DiscreteFunction<SpaceType, VectorType> DiscreteFunctionType;
    VectorType source_vector(space.mapper().size());
    DiscreteFunctionType source(space, source_vector);
    for (auto entity_ptr = grid_part_view.template begin<0>(); entity_ptr != grid_part_view.template end<0>();
         ++entity_ptr) {
      const auto& entity = *entity_ptr;
      const auto center  = entity.geometry().center();
      double value = 0.0;
      if (center[0] > 0.5)
        value                      = 1.0;
      auto local_source            = source.local_discrete_function(entity);
      auto local_source_DoF_vector = local_source->vector();
      for (size_t local_DoF = 0; local_DoF < local_source_DoF_vector.size(); ++local_DoF)
        local_source_DoF_vector.set(local_DoF, value);
    }
    // apply operator
    VectorType range_vector(space.mapper().size());
    DiscreteFunctionType range(space, range_vector);
    Operators::OswaldInterpolation<typename SpaceType::GridViewType> oswald_operator(space.grid_view());
    oswald_operator.apply(source, range);
    // TODO: test result
  } // ... produces_correct_results()
}; // struct Oswald_Interpolation_Operator


#endif // DUNE_GDT_TEST_OPERATORS_OSWALDINTERPOLATION_HH
