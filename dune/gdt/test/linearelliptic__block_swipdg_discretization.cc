#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/projections.hh>

#include <dune/gdt/test/linearelliptic/discretizers/block-ipdg.hh>
#include <dune/gdt/test/linearelliptic/problems/ESV2007.hh>

using namespace Dune;

GTEST_TEST(linearelliptic_block_SWIPDG_discretization, coincides_with_SWIPDG)
{
  typedef Dune::ALUGrid<2, 2, simplex, conforming> GridType;
  size_t inner_boundary_index = std::numeric_limits<size_t>::max() - 42;
  auto grid_provider = XT::Grid::make_cube_dd_subdomains_grid<GridType>(
      {-1, -1}, {1, 1}, {9, 9}, 1, {0, 0}, {3, 3}, 0, inner_boundary_index);

  XT::Common::Configuration local_boundary_cfg;
  local_boundary_cfg["type"] = "xt.grid.boundaryinfo.boundarysegment";
  local_boundary_cfg["default"] = "dirichlet";
  local_boundary_cfg["neumann"] =
      "[" + XT::Common::to_string(inner_boundary_index) + " " + XT::Common::to_string(inner_boundary_index + 1) + "]";

  grid_provider.visualize_dd("grid", local_boundary_cfg);

  GDT::LinearElliptic::ESV2007TestCase<GridType> test_case;
  const auto& problem = test_case.problem();

  auto block_ipdg_disc =
      GDT::LinearElliptic::BlockIpdgDiscretizer<GridType>::discretize(grid_provider, problem, -1, inner_boundary_index);
  auto block_ipdg_solution_vector = block_ipdg_disc.create_vector();
  block_ipdg_disc.solve(block_ipdg_solution_vector);

  auto ipdg_disc =
      GDT::LinearElliptic::IpdgDiscretizer<GridType, XT::Grid::Layers::leaf, GDT::ChooseSpaceBackend::fem>::discretize(
          grid_provider, problem);
  auto ipdg_solution_vector = ipdg_disc.create_vector();
  ipdg_disc.solve(ipdg_solution_vector);
  ipdg_disc.visualize(ipdg_solution_vector, "ipdg_solution", "solution");

  auto block_ipdg_solution_in_block_ipdg_space =
      GDT::make_const_discrete_function(block_ipdg_disc.ansatz_space(), block_ipdg_solution_vector);
  auto block_ipdg_solution_in_ipdg_space =
      GDT::make_discrete_function<decltype(block_ipdg_solution_vector)>(ipdg_disc.ansatz_space(), "solution");

  GDT::project(block_ipdg_solution_in_block_ipdg_space, block_ipdg_solution_in_ipdg_space);
  block_ipdg_solution_in_ipdg_space.visualize("block_ipdg_solution");

  block_ipdg_solution_in_ipdg_space.vector() -= ipdg_solution_vector;
  for (auto& element : block_ipdg_solution_in_ipdg_space.vector())
    element = std::abs(element);
  block_ipdg_solution_in_ipdg_space.visualize("difference");

  EXPECT_LE(block_ipdg_solution_in_ipdg_space.vector().sup_norm(), 1e-15);
}
