// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
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

#include <boost/filesystem.hpp>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/misc/gridwidth.hh>
#include <dune/fem/gridpart/levelgridpart.hh>

#include <dune/stuff/common/parameter/tree.hh>
#include <dune/stuff/common/logging.hh>
#include <dune/stuff/grid/provider/cube.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/functions/expression.hh>

#include "discretizations.hh"

using namespace Dune;
using namespace Dune::GDT;


int main(int argc, char** argv)
{
  try {
    // defines
    static const unsigned int dimDomain = 2;
    static const unsigned int dimRange  = 1;
    typedef double RangeFieldType;
    const std::string static_id = "example.elliptic.convergence_study";

    // read or write settings file
    if (!boost::filesystem::exists(static_id + ".settings")) {
      std::cout << "writing default settings to '" << static_id + ".settings'" << std::endl;
      std::ofstream file;
      file.open(static_id + ".settings");
      file << "integration_order = 4" << std::endl;
      file << "num_refinements   = 3" << std::endl;
      file.close();
    } else {
      std::cout << "reading default settings from'" << static_id + ".settings'" << std::endl;
      std::cout << std::endl;
      Stuff::Common::ExtendedParameterTree settings(argc, argv, static_id + ".settings");
      const size_t integration_order = settings.get("integration_order", 4);
      const size_t num_refinements = settings.get("num_refinements", 3);
      DSC_LOG.create(Stuff::Common::LOG_INFO, "", "");
      auto& info = DSC_LOG.info();

      info << "==================================================================================" << std::endl;
      info << "  convergence study (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)" << std::endl;
      info << "  dimDomain = " << dimDomain << ", dimRange = " << dimRange << std::endl;
      info << "  integration_order = " << integration_order << std::endl;
      info << "  num_refinements   = " << num_refinements << std::endl;
      info << "==================================================================================" << std::endl;

      // mpi
      Dune::Fem::MPIManager::initialize(argc, argv);

      // prepare the grid
      typedef Dune::ALUConformGrid<dimDomain, dimDomain> GridType;
      typedef typename GridType::ctype DomainFieldType;
      typedef Stuff::GridProviderCube<GridType> GridProviderType;
      GridProviderType grid_provider(-1, 1, 2);
      GridType& grid = *(grid_provider.grid());
      grid.globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 2); ++ii)
        grid.globalRefine(GridType::refineStepsForHalf);

      // prepare the analytical functions
      typedef Stuff::FunctionConstant<DomainFieldType, dimDomain, RangeFieldType, dimRange> ConstantFunctionType;
      typedef Stuff::FunctionExpression<DomainFieldType, dimDomain, RangeFieldType, dimRange> ExpressionFunctionType;
      const ConstantFunctionType diffusion(1.0);
      const ExpressionFunctionType force(
          "x", "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])", integration_order);
      const ConstantFunctionType dirichlet(0.0);
      const ExpressionFunctionType exact_solution("x",
                                                  "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                                                  {"-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
                                                   "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"},
                                                  integration_order);

      // compute exact solution norm
      typedef typename Fem::LevelGridPart<GridType> GridPartType;
      const GridPartType finest_grid_part(grid, 2 * num_refinements + 3);
      const Stuff::GridboundaryAllDirichlet<typename GridPartType::GridViewType> boundary_info;
      ProductOperator::L2<GridPartType> l2_product_operator(finest_grid_part);
      const RangeFieldType exact_solution_l2_norm =
          std::sqrt(l2_product_operator.apply2(exact_solution, exact_solution));

      // prepare the data structures
      std::vector<size_t> num_grid_elements(num_refinements + 1);
      std::vector<double> grid_width(num_refinements + 1);
      std::vector<RangeFieldType> cg_l2_absolute_errors(num_refinements + 1);
      std::vector<RangeFieldType> cg_l2_relative_errors(num_refinements + 1);
      // print table header
      info << "continuous galerkin, polOrder = 1, error against exact solution" << std::endl;
      info << "=================================================================" << std::endl;
      info << "        grid         |   L2 (absolute)     |    L2 (relative)" << std::endl;
      info << "---------------------+---------------------+---------------------" << std::endl;
      info << "     size |    width |    error |      EOC |    error |      EOC" << std::endl;
      info << "==========+==========+==========+==========+==========+==========" << std::endl;

      // iterate
      for (size_t ii = 0; ii <= num_refinements; ++ii) {
        if (ii > 0)
          info << "----------+----------+----------+----------+----------+----------" << std::endl;

        const size_t level = 2 * ii + 1;
        GridPartType level_grid_part(grid, level);
        num_grid_elements[ii] = grid.size(level, 0);
        grid_width[ii] = Fem::GridWidth::calcGridWidth(level_grid_part);
        info << " " << std::setw(8) << num_grid_elements[ii] << " | " << std::setw(8) << std::setprecision(2)
             << std::scientific << grid_width[ii] << " | " << std::flush;

        // discretize
        Example::CGDiscretization<GridPartType, 1> cg_discretization(
            level_grid_part, boundary_info, diffusion, force, dirichlet);
        auto cg_solution = cg_discretization.solve();
        cg_discretization.visualize(*(cg_solution->vector()),
                                    static_id + ".cg_solution." + Stuff::Common::toString(ii));

        // compute errors
        auto cg_errors            = cg_discretization.compute_errors(exact_solution, *(cg_solution));
        cg_l2_absolute_errors[ii] = cg_errors[0];
        info << std::setw(8) << std::setprecision(2) << std::scientific << cg_l2_absolute_errors[ii] << " | "
             << std::flush;
        if (ii == 0)
          info << std::setw(11) << "---- | " << std::flush;
        else
          info << std::setw(8) << std::setprecision(2) << std::fixed
               << std::log(cg_l2_absolute_errors[ii] / cg_l2_absolute_errors[ii - 1])
                      / std::log(grid_width[ii] / grid_width[ii - 1])
               << " | " << std::flush;

        cg_l2_relative_errors[ii] = cg_l2_absolute_errors[ii] / exact_solution_l2_norm;
        info << std::setw(8) << std::setprecision(2) << std::scientific << cg_l2_relative_errors[ii] << " | "
             << std::flush;
        if (ii == 0)
          info << std::setw(8) << "----" << std::flush;
        else
          info << std::setw(8) << std::setprecision(2) << std::fixed
               << std::log(cg_l2_relative_errors[ii] / cg_l2_relative_errors[ii - 1])
                      / std::log(grid_width[ii] / grid_width[ii - 1])
               << std::flush;

        info << std::endl;

      } // iterate
    } // read or write settings file

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
