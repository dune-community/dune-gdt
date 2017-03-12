#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/test/linearelliptic/discretizers/block-ipdg.hh>
#include <dune/gdt/test/linearelliptic/problems/ESV2007.hh>

using namespace Dune;

GTEST_TEST(linearelliptic_block_SWIPDG_discretization, coincides_with_SWIPDG)
{
  typedef Dune::ALUGrid<2, 2, simplex, conforming> GridType;
  size_t inner_boundary_index = std::numeric_limits<size_t>::max() - 42;
  auto grid_provider = XT::Grid::make_cube_dd_subdomains_grid<GridType>(
      {0, 0}, {1, 1}, {9, 9}, 1, {0, 0}, {3, 3}, 0, inner_boundary_index);
  grid_provider.visualize_dd("grid", true);

  GDT::LinearElliptic::ESV2007TestCase<GridType> test_case;
  const auto& problem = test_case.problem();

  auto block_ipdg_disc =
      GDT::LinearElliptic::BlockIpdgDiscretizer<GridType>::discretize(grid_provider, problem, -1, inner_boundary_index);
  auto block_ipdg_solution = block_ipdg_disc.create_vector();
  block_ipdg_disc.solve(block_ipdg_solution);

  auto ipdg_disc =
      GDT::LinearElliptic::IpdgDiscretizer<GridType, XT::Grid::Layers::leaf, GDT::ChooseSpaceBackend::fem>::discretize(
          grid_provider, problem);
  auto ipdg_solution = ipdg_disc.create_vector();
  ipdg_disc.solve(ipdg_solution);

  EXPECT_LE((ipdg_solution - block_ipdg_solution).sup_norm(), 1e-15);
}
