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
#include <dune/stuff/functions/checkerboard.hh>
#include <dune/stuff/functions/visualization.hh>

#include <dune/gdt/operator/prolongations.hh>

#include "discretizations.hh"

using namespace Dune;
using namespace Dune::GDT;


template <class DiscretizationType>
class ConvergenceStudy
{
public:
  template <class GridType, class ReferenceDiscretizationType, class ReferenceSolutionType>
  static void run(GridType& grid, const ReferenceDiscretizationType& reference_discretization,
                  const size_t num_refinements, const ReferenceSolutionType& reference_solution,
                  const std::string header, const std::string plot_prefix)
  {
    typedef typename ReferenceDiscretizationType::DomainFieldType DomainFieldType;
    typedef typename ReferenceDiscretizationType::RangeFieldType RangeFieldType;
    typedef typename ReferenceDiscretizationType::SpaceType::GridPartType GridPartType;
    // prepare the data structures
    std::vector<size_t> num_grid_elements(num_refinements + 1);
    std::vector<DomainFieldType> grid_width(num_refinements + 1);
    std::vector<RangeFieldType> l2_errors(num_refinements + 1);
    std::vector<RangeFieldType> h1_errors(num_refinements + 1);
    // print table header
    const std::string bar = "=================================================================";
    auto& info            = DSC_LOG.info();
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
    ProductOperator::L2<GridPartType> l2_product_operator(reference_discretization.grid_part());
    const RangeFieldType reference_solution_l2_norm =
        std::sqrt(l2_product_operator.apply2(reference_solution, reference_solution));
    ProductOperator::H1<GridPartType> h1_product_operator(reference_discretization.grid_part());
    const RangeFieldType reference_solution_h1_error =
        std::sqrt(h1_product_operator.apply2(reference_solution, reference_solution));

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
      const DiscretizationType discretization(level_grid_part,
                                              reference_discretization.boundary_info(),
                                              reference_discretization.diffusion(),
                                              reference_discretization.force(),
                                              reference_discretization.dirichlet());
      auto discrete_solution = discretization.solve();
      discretization.visualize(*(discrete_solution->vector()),
                               plot_prefix + "." + discretization.id() + ".solution." + Stuff::Common::toString(ii));
      typedef typename ReferenceDiscretizationType::DiscreteFunctionType DiscreteFunctionType;
      typedef typename ReferenceDiscretizationType::VectorType VectorType;
      auto discrete_solution_on_reference_grid = std::make_shared<DiscreteFunctionType>(
          reference_discretization.space(),
          std::make_shared<VectorType>(reference_discretization.space().mapper().size()),
          "discrete_solution");
      ProlongationOperator::Generic::apply(*discrete_solution, *discrete_solution_on_reference_grid);

      // compute errors
      const auto errors =
          reference_discretization.compute_errors(reference_solution, *(discrete_solution_on_reference_grid));
      // * L2
      l2_errors[ii] = errors[0];
      info << std::setw(8) << std::setprecision(2) << std::scientific << l2_errors[ii] / reference_solution_l2_norm
           << " | " << std::flush;
      if (ii == 0)
        info << std::setw(11) << "---- | " << std::flush;
      else
        info << std::setw(8) << std::setprecision(2) << std::fixed
             << std::log(l2_errors[ii] / l2_errors[ii - 1]) / std::log(grid_width[ii] / grid_width[ii - 1]) << " | "
             << std::flush;
      // * H1
      h1_errors[ii] = errors[1];
      info << std::setw(8) << std::setprecision(2) << std::scientific << h1_errors[ii] / reference_solution_h1_error
           << " | " << std::flush;
      if (ii == 0)
        info << std::setw(8) << "----" << std::flush;
      else
        info << std::setw(8) << std::setprecision(2) << std::fixed
             << std::log(h1_errors[ii] / h1_errors[ii - 1]) / std::log(grid_width[ii] / grid_width[ii - 1])
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
      std::cout << "reading default settings from'" << static_id + ".settings':" << std::endl;
      Stuff::Common::ExtendedParameterTree settings(argc, argv, static_id + ".settings");
      const size_t integration_order = settings.get("integration_order", 4);
      const size_t num_refinements = settings.get("num_refinements", 3);
      DSC_LOG.create(Stuff::Common::LOG_INFO, "", "");
      auto& info = DSC_LOG.info();
      info << "  dimDomain = " << dimDomain << ", dimRange = " << dimRange << std::endl;
      info << "  integration_order = " << integration_order << std::endl;
      info << "  num_refinements   = " << num_refinements << std::endl;
      info << std::endl;

      // some preparations
      typedef Dune::ALUConformGrid<dimDomain, dimDomain> GridType;
      typedef Stuff::GridProviderCube<GridType> GridProviderType;
      typedef typename Fem::LevelGridPart<GridType> GridPartType;
      typedef VTKWriter<typename GridPartType::GridViewType> VTKWriterType;
      typedef typename GridType::ctype DomainFieldType;
      typedef FieldVector<DomainFieldType, dimDomain> DomainType;
      typedef Stuff::FunctionConstant<DomainFieldType, dimDomain, RangeFieldType, dimRange> ConstantFunctionType;
      typedef Stuff::FunctionExpression<DomainFieldType, dimDomain, RangeFieldType, dimRange> ExpressionFunctionType;
      typedef Stuff::FunctionCheckerboard<DomainFieldType, dimDomain, RangeFieldType, dimRange>
          CheckerboardFunctionType;
      typedef Stuff::FunctionVisualization<typename GridPartType::GridViewType,
                                           DomainFieldType,
                                           dimDomain,
                                           RangeFieldType,
                                           dimRange> VisualizationFunctionType;
      typedef Example::CGDiscretization<GridPartType, 1> CG_1_DiscretizationType;
      typedef Example::SIPDGDiscretization<GridPartType, 1> SIPDG_1_DiscretizationType;
      Dune::Fem::MPIManager::initialize(argc, argv);

      info << "============================================================================" << std::endl;
      info << "  Testcase 1: smooth data, homogeneous dirichlet" << std::endl;
      info << "              (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)" << std::endl;
      info << "  domain = [-1, 1] x [-1 , 1]" << std::endl;
      info << "  diffusion = 1" << std::endl;
      info << "  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)" << std::endl;
      info << "  dirichlet = 0" << std::endl;
      info << "  exact solution = cos(1/2 pi x) cos(1/2 pi y)" << std::endl;
      info << "============================================================================" << std::endl;
      auto grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(-1, 1, 2));
      auto grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 2); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      auto reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      auto vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      const Stuff::GridboundaryAllDirichlet<typename GridPartType::GridViewType> testcase_1_boundary_info;
      const ConstantFunctionType testcase_1_diffusion(1.0);
      const ExpressionFunctionType testcase_1_force(
          "x", "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])", integration_order);
      const ConstantFunctionType testcase_1_dirichlet(0.0);
      const ExpressionFunctionType testcase_1_exact_solution(
          "x",
          "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
          {"-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
           "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"},
          integration_order);
      vtk_writer->addCellData(std::make_shared<VisualizationFunctionType>(testcase_1_force));
      vtk_writer->write("testcase_1.force");
      vtk_writer->clear();
      vtk_writer->addCellData(std::make_shared<VisualizationFunctionType>(testcase_1_exact_solution));
      vtk_writer->write("testcase_1.exact_solution");
      vtk_writer->clear();

      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_1_cg_1_reference_discretization(
          *reference_grid_part, testcase_1_boundary_info, testcase_1_diffusion, testcase_1_force, testcase_1_dirichlet);
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_1_cg_1_reference_discretization,
                                                     num_refinements,
                                                     testcase_1_exact_solution,
                                                     "continuous galerkin, polOrder 1",
                                                     "testcase_1");
      info << std::endl;
      // symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_1_DiscretizationType testcase_1_sipdg_1_reference_discretization(
          *reference_grid_part, testcase_1_boundary_info, testcase_1_diffusion, testcase_1_force, testcase_1_dirichlet);
      ConvergenceStudy<SIPDG_1_DiscretizationType>::run(*grid,
                                                        testcase_1_sipdg_1_reference_discretization,
                                                        num_refinements,
                                                        testcase_1_exact_solution,
                                                        "symmetric interior penalty discontinuous galerkin, polOrder 1",
                                                        "testcase_1");
      info << std::endl;

      info << "=========================================================" << std::endl;
      info << "  Testcase 2: smooth data, smooth dirichlet" << std::endl;
      info << "              (see page 858 in Epshteyn, Riviere, 2007)" << std::endl;
      info << "  domain = [0, 1] x [0 , 1]" << std::endl;
      info << "  diffusion = 1" << std::endl;
      info << "  force     = 64 pi^2 (cos(8  pi x) + cos(8 pi y))" << std::endl;
      info << "  dirichlet = cos(8 pi x) + cos(8 pi y)" << std::endl;
      info << "  exact solution = cos(8 pi x) + cos(8 pi y)" << std::endl;
      info << "=========================================================" << std::endl;
      grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(0, 1, 8));
      grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 2); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      const Stuff::GridboundaryAllDirichlet<typename GridPartType::GridViewType> testcase_2_boundary_info;
      const ConstantFunctionType testcase_2_diffusion(1.0);
      const ExpressionFunctionType testcase_2_force(
          "x", "64.0 * pi * pi * (cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1]))", integration_order);
      const ExpressionFunctionType testcase_2_dirichlet(
          "x", "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])", integration_order);
      const ExpressionFunctionType testcase_2_exact_solution(
          "x",
          "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])",
          {"-8.0 * pi * sin(8.0 * pi * x[0])", "-8.0 * pi * sin(8.0 * pi * x[1])"},
          integration_order);
      vtk_writer->addCellData(std::make_shared<VisualizationFunctionType>(testcase_2_force));
      vtk_writer->write("testcase_2.force");
      vtk_writer->clear();
      vtk_writer->addCellData(std::make_shared<VisualizationFunctionType>(testcase_2_exact_solution));
      vtk_writer->write("testcase_2.exact_solution");
      vtk_writer->clear();

      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_2_cg_1_reference_discretization(
          *reference_grid_part, testcase_2_boundary_info, testcase_2_diffusion, testcase_2_force, testcase_2_dirichlet);
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_2_cg_1_reference_discretization,
                                                     num_refinements,
                                                     testcase_2_exact_solution,
                                                     "continuous galerkin, polOrder 1",
                                                     "testcase_2");
      info << std::endl;
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
