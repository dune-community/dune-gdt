// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#ifdef HAVE_ALUGRID_SERIAL_H
  #define ENABLE_ALUGRID 1
  #undef HAVE_ALUGRID
  #define HAVE_ALUGRID 1
  #include <dune/grid/alugrid.hh>
#else
  static_assert(false, "This study requires a serial alugrid!");
#endif // HAVE_ALUGRID_SERIAL_H

#include <memory>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/gridpart/levelgridpart.hh>

#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include "discretizations.hh"

using namespace Dune;


int main(int argc, char** argv)
{
  try {
    // defines
    static const unsigned int DUNE_UNUSED(dimDomain) = 2;
    static const unsigned int DUNE_UNUSED(dimRange) = 1;
    typedef double RangeFieldType;
    const size_t num_refinements = 5;
    const size_t integration_order = 4;
    const std::string static_id = "example.elliptic.convergence_study";

    // mpi
    Dune::Fem::MPIManager::initialize(argc, argv);

    // prepare the grid
    typedef Dune::ALUConformGrid< dimDomain, dimDomain > GridType;
    typedef typename GridType::ctype DomainFieldType;
    typedef Stuff::GridProviderCube< GridType > GridProviderType;
    GridProviderType grid_provider(-1, 1, 2);
    GridType& grid = *(grid_provider.grid());
    grid.globalRefine(1);
    grid_provider.visualize(static_id + ".0");
    for (size_t ii = 1; ii <= num_refinements; ++ii) {
      grid.globalRefine(GridType::refineStepsForHalf);
      grid_provider.visualize(static_id + "." + Stuff::Common::toString(ii));
    }

    // prepare the analytical functions
    typedef Stuff::FunctionConstant< DomainFieldType, dimDomain, RangeFieldType, dimRange >   ConstantFunctionType;
    typedef Stuff::FunctionExpression< DomainFieldType, dimDomain, RangeFieldType, dimRange > ExpressionFunctionType;
    const ConstantFunctionType    diffusion(1.0);
    const ExpressionFunctionType  force("x",
                                        "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                                        integration_order);
    const ConstantFunctionType    dirichlet(0.0);
    const ExpressionFunctionType  exact_solution("x",
                                                 "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                                                 {"-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                                                  "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"},
                                                 integration_order);

    // discretize
    typedef typename Fem::LevelGridPart< GridType > GridPartType;
    GridPartType level_grid_part(grid, 1);
    const Stuff::GridboundaryAllDirichlet< typename GridPartType::GridViewType > boundary_info;
    Example::CGDiscretization< GridPartType, 1 > cg_discretization(level_grid_part,
                                                                   boundary_info,
                                                                   diffusion, force, dirichlet);
    auto cg_solution = cg_discretization.solve();
    cg_discretization.visualize(*(cg_solution->vector()), static_id + ".cg_solution.1");
    cg_discretization.compute_errors(exact_solution, *(cg_solution));

    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
  } // try
} // main
