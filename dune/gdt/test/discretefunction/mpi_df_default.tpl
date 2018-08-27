// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2017 - 2018)

#include <dune/xt/common/test/main.hxx> // <- This one has to come first!

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/gridprovider.hh>

#include <dune/gdt/test/grids.hh>
#include <dune/gdt/test/prolongations/prolongations.hh>
#include <dune/gdt/test/prolongations/l2.hh>
#include <dune/gdt/spaces/cg.hh>

#include <dune/gdt/spaces/dg.hh>
#include <dune/gdt/spaces/fv.hh>
#include <dune/gdt/spaces/rt/default.hh>


using namespace Dune::GDT::Test;
using namespace Dune::GDT;

// clang-format off
{% for SpaceType,Name,View in config.spaces_with_names %}


GTEST_TEST(DiscretefunctionTest_{{Name}}, visualize)
{
  using Grid = Dune::XT::Grid::extract_grid_t<{{View}}>;
  auto grid = Dune::XT::Grid::make_cube_grid<Grid>(0.0, 1.0, 6u);
  grid.grid().globalRefine(1);
  boost::format viz_name{"df_viz_%s_lvl_%d"};

  for (auto&& ii : Dune::XT::Common::value_range(grid.max_level()+1)) {
    {{View}} view = grid.level_view(ii);
    {{SpaceType}} space{view};
    DiscreteFunction<{{SpaceType}}> function{space};
    const std::string view_name{"{{Name}}"};
    function.visualize((viz_name % view_name % ii).str());
    auto func_copy = function;
    func_copy.visualize((viz_name % (view_name + "_copy") % ii).str());
  }
}


{% endfor %}
// clang-format on