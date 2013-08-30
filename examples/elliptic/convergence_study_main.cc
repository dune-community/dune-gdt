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

#include <dune/gdt/operator/prolongations.hh>

#include "discretizations.hh"

using namespace Dune;
using namespace Dune::GDT;


template< class DiscretizationType >
class ConvergenceStudy
{
public:
  template< class GridType, class ReferenceDiscretizationType, class ReferenceSolutionType >
  static void run(GridType& grid,
                  const ReferenceDiscretizationType& reference_discretization,
                  const size_t num_refinements,
                  const ReferenceSolutionType& reference_solution,
                  const std::string header)
  {
    typedef typename ReferenceDiscretizationType::DomainFieldType         DomainFieldType;
    typedef typename ReferenceDiscretizationType::RangeFieldType          RangeFieldType;
    typedef typename ReferenceDiscretizationType::SpaceType::GridPartType GridPartType;
    // prepare the data structures
    std::vector< size_t > num_grid_elements(num_refinements + 1);
    std::vector< DomainFieldType > grid_width(num_refinements + 1);
    std::vector< RangeFieldType > l2_errors(num_refinements + 1);
    std::vector< RangeFieldType > h1_errors(num_refinements + 1);
    // print table header
    const std::string bar = "=================================================================";
    auto& info = DSC_LOG.info();
    info << header << std::endl;
    if (header.size() > bar.size())
      info << Stuff::Common::whitespaceify(header, '=') << std::endl;
    else
      info << bar << std::endl;
    info << "        grid         |    L2 (relative)    |    H1 (relative)    " << std::endl;
    info << "---------------------+---------------------+---------------------" << std::endl;
    info << "     size |    width |    error |      EOC |    error |      EOC " << std::endl;
    info << "==========+==========+==========+==========+==========+==========" << std::endl;

    // compute norm of reference solution
    ProductOperator::L2< GridPartType > l2_product_operator(reference_discretization.grid_part());
    const RangeFieldType reference_solution_l2_norm = std::sqrt(l2_product_operator.apply2(reference_solution,
                                                                                           reference_solution));
    ProductOperator::H1< GridPartType > h1_product_operator(reference_discretization.grid_part());
    const RangeFieldType reference_solution_h1_error = std::sqrt(h1_product_operator.apply2(reference_solution,
                                                                                            reference_solution));

    // iterate
    for (size_t ii = 0; ii <= num_refinements; ++ii) {
      if (ii > 0)
        info << "----------+----------+----------+----------+----------+----------" << std::endl;

      const size_t level = 2*ii + 1;
      GridPartType level_grid_part(grid, level);
      num_grid_elements[ii] = grid.size(level, 0);
      grid_width[ii] = Fem::GridWidth::calcGridWidth(level_grid_part);
      info << " " << std::setw(8) << num_grid_elements[ii]
              << " | " << std::setw(8) << std::setprecision(2) << std::scientific << grid_width[ii]
              << " | " << std::flush;

      // discretize
      const DiscretizationType discretization(level_grid_part,
                                              reference_discretization.boundary_info(),
                                              reference_discretization.diffusion(),
                                              reference_discretization.force(),
                                              reference_discretization.dirichlet());
      auto discrete_solution = discretization.solve();
      typedef typename ReferenceDiscretizationType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename ReferenceDiscretizationType::VectorType VectorType;
      auto discrete_solution_on_reference_grid
          = std::make_shared< DiscreteFunctionType >(reference_discretization.space(),
                                                     std::make_shared< VectorType >(reference_discretization.space().mapper().size()),
                                                     "discrete_solution");
      ProlongationOperator::Generic::apply(*discrete_solution, *discrete_solution_on_reference_grid);

      // compute errors
      const auto errors = reference_discretization.compute_errors(reference_solution,
                                                                  *(discrete_solution_on_reference_grid));
      // * L2
      l2_errors[ii] = errors[0];
      info << std::setw(8) << std::setprecision(2) << std::scientific << l2_errors[ii] / reference_solution_l2_norm
              << " | " << std::flush;
      if (ii == 0)
        info << std::setw(11) << "---- | " << std::flush;
      else
        info << std::setw(8) << std::setprecision(2) << std::fixed
             << std::log(l2_errors[ii] / l2_errors[ii - 1])
                / std::log(grid_width[ii] / grid_width[ii - 1])
            << " | " << std::flush;
      // * H1
      h1_errors[ii] = errors[1];
      info << std::setw(8) << std::setprecision(2) << std::scientific << h1_errors[ii] / reference_solution_h1_error
              << " | " << std::flush;
      if (ii == 0)
        info << std::setw(8) << "----" << std::flush;
      else
        info << std::setw(8) << std::setprecision(2) << std::fixed
             << std::log(h1_errors[ii] / h1_errors[ii - 1])
                / std::log(grid_width[ii] / grid_width[ii - 1])
            << std::flush;

      info << std::endl;
    } // iterate
  }
}; // class ConvergenceStudy


int main(int argc, char** argv)
{
  try {
    // defines
    static const unsigned int dimDomain = 2;
    static const unsigned int dimRange = 1;
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
      const size_t num_refinements   = settings.get("num_refinements",   3);
      DSC_LOG.create(Stuff::Common::LOG_INFO, "", "");
      auto& info = DSC_LOG.info();

      info << "==================================================================================" << std::endl;
      info << "  convergence study (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)"   << std::endl;
      info << "  dimDomain = " << dimDomain << ", dimRange = " << dimRange << std::endl;
      info << "  integration_order = " << integration_order << std::endl;
      info << "  num_refinements   = " << num_refinements << std::endl;
      info << "==================================================================================" << std::endl;

      // mpi
      Dune::Fem::MPIManager::initialize(argc, argv);

      // prepare the grid
      typedef Dune::ALUConformGrid< dimDomain, dimDomain > GridType;
      typedef typename GridType::ctype DomainFieldType;
      typedef Stuff::GridProviderCube< GridType > GridProviderType;
      GridProviderType grid_provider(-1, 1, 2);
      GridType& grid = *(grid_provider.grid());
      grid.globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 2); ++ii)
        grid.globalRefine(GridType::refineStepsForHalf);

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
      typedef typename Fem::LevelGridPart< GridType > GridPartType;
      const Stuff::GridboundaryAllDirichlet< typename GridPartType::GridViewType > boundary_info;

      // compute reference
      const GridPartType reference_grid_part(grid, 2*num_refinements + 3);

      // continuous galerkin discretization
      typedef Example::CGDiscretization< GridPartType, 1 > CG_1_DiscretizationType;
      const CG_1_DiscretizationType cg_1_reference_discretization(reference_grid_part, boundary_info, diffusion, force, dirichlet);
      ConvergenceStudy< CG_1_DiscretizationType >::run(grid,
                                                       cg_1_reference_discretization,
                                                       num_refinements,
                                                       exact_solution,
                                                       "continuous galerkin, polOrder = 1, error against exact solution");
      info << std::endl;
      typedef Example::CGDiscretization< GridPartType, 2 > CG_2_DiscretizationType;
      const CG_2_DiscretizationType cg_2_reference_discretization(reference_grid_part, boundary_info, diffusion, force, dirichlet);
      ConvergenceStudy< CG_2_DiscretizationType >::run(grid,
                                                       cg_2_reference_discretization,
                                                       num_refinements,
                                                       exact_solution,
                                                       "continuous galerkin, polOrder = 2, error against exact solution");
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
