// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if !defined(HAVE_ALUGRID)
static_assert(false, "This study requires alugrid!");
#endif

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
#include <dune/stuff/functions/spe10.hh>

#include <dune/gdt/operator/prolongations.hh>

#include "discretization-cg.hh"
#include "discretization-sipdg.hh"

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
                                              reference_discretization.dirichlet(),
                                              reference_discretization.neumann());
      try {
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
      } catch (Dune::MathError&) {
        info << Stuff::Common::colorStringRed("ERROR:") << " linear solver failed!" << std::endl;
      }
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
      file << "spe_10_model_1_datafile = perm_case1.dat" << std::endl;
      file.close();
    } else {
      std::cout << "reading default settings from '" << static_id + ".settings':" << std::endl;
      Stuff::Common::ExtendedParameterTree settings(argc, argv, static_id + ".settings");
      const size_t integration_order            = settings.get("integration_order", 4);
      const size_t num_refinements              = settings.get("num_refinements", 3);
      const std::string spe_10_model_1_datafile = settings.get("spe_10_model_1_datafile", "perm_case1.dat");
      DSC_LOG.create(Stuff::Common::LOG_INFO, "", "", "");
      auto& info = DSC_LOG.info();
      info << "  dimDomain = " << dimDomain << ", dimRange = " << dimRange << std::endl;
      info << "  integration_order = " << integration_order << std::endl;
      info << "  num_refinements   = " << num_refinements << std::endl;
      const bool spe_10_model_1_datafile_found = boost::filesystem::exists(spe_10_model_1_datafile);
      if (spe_10_model_1_datafile_found)
        info << "  spe_10_model_1_datafile: " << spe_10_model_1_datafile << std::endl;
      else
        info << Stuff::Common::colorString("WARNING:") << " spe_10_model_1_datafile '" << spe_10_model_1_datafile
             << "' "
             << "could not be found!" << std::endl;
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
      typedef Stuff::FunctionSpe10Model1<DomainFieldType, dimDomain, RangeFieldType, dimRange> SPE10Model1FunctionType;
      typedef Discretization::ContinuousGalerkin<GridPartType, 1> CG_1_DiscretizationType;
      typedef Discretization::SymmetricInteriorPenaltyDG<GridPartType, 1> SIPDG_1_DiscretizationType;
      typedef Discretization::SymmetricInteriorPenaltyDGAdhoc<GridPartType, 1> SIPDG_Adhoc_1_DiscretizationType;
      Dune::Fem::MPIManager::initialize(argc, argv);

      info << "+==============================================================================+" << std::endl;
      info << "|+============================================================================+|" << std::endl;
      info << "||  Testcase 1: smooth data, homogeneous dirichlet                            ||" << std::endl;
      info << "||              (see testcase 1, page 23 in Ern, Stephansen, Vohralik, 2007)  ||" << std::endl;
      info << "|+----------------------------------------------------------------------------+|" << std::endl;
      info << "||  domain = [-1, 1] x [-1 , 1]                                               ||" << std::endl;
      info << "||  diffusion = 1                                                             ||" << std::endl;
      info << "||  force     = 1/2 pi^2 cos(1/2 pi x) cos(1/2 pi y)                          ||" << std::endl;
      info << "||  dirichlet = 0                                                             ||" << std::endl;
      info << "||  exact solution = cos(1/2 pi x) cos(1/2 pi y)                              ||" << std::endl;
      info << "|+============================================================================+|" << std::endl;
      info << "+==============================================================================+" << std::endl;
      info << std::endl;
      auto grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(-1, 1, 2));
      auto grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 1); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      auto reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      auto vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      const Stuff::GridboundaryAllDirichlet<typename GridPartType::IntersectionType> testcase_1_boundary_info;
      const ConstantFunctionType testcase_1_diffusion(1.0);
      const ExpressionFunctionType testcase_1_force(
          "x", "0.5 * pi * pi * cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])", integration_order);
      const ConstantFunctionType testcase_1_dirichlet(0.0);
      const ConstantFunctionType testcase_1_neumann(0.0);
      const ExpressionFunctionType testcase_1_exact_solution(
          "x",
          "cos(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
          {"-0.5 * pi * sin(0.5 * pi * x[0]) * cos(0.5 * pi * x[1])",
           "-0.5 * pi * cos(0.5 * pi * x[0]) * sin(0.5 * pi * x[1])"},
          integration_order);

      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_1_cg_1_reference_discretization(*reference_grid_part,
                                                                             testcase_1_boundary_info,
                                                                             testcase_1_diffusion,
                                                                             testcase_1_force,
                                                                             testcase_1_dirichlet,
                                                                             testcase_1_neumann);
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_1_cg_1_reference_discretization,
                                                     num_refinements,
                                                     testcase_1_exact_solution,
                                                     "continuous galerkin discretization, polOrder 1",
                                                     "testcase_1");
      info << std::endl;
      // symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_1_DiscretizationType testcase_1_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                   testcase_1_boundary_info,
                                                                                   testcase_1_diffusion,
                                                                                   testcase_1_force,
                                                                                   testcase_1_dirichlet,
                                                                                   testcase_1_neumann);
      ConvergenceStudy<SIPDG_1_DiscretizationType>::run(*grid,
                                                        testcase_1_sipdg_1_reference_discretization,
                                                        num_refinements,
                                                        testcase_1_exact_solution,
                                                        "SIP discontinuous galerkin discretization, polOrder 1",
                                                        "testcase_1");
      info << std::endl;
      //      const Example::NewSIPDGDiscretization< GridPartType, 2 >
      //          testcase_1_new_sipdg_2_reference_discretization(*reference_grid_part,
      //                                                          testcase_1_boundary_info,
      //                                                          testcase_1_diffusion,
      //                                                          testcase_1_force,
      //                                                          testcase_1_dirichlet,
      //                                                          testcase_1_neumann);
      //      ConvergenceStudy< Example::NewSIPDGDiscretization< GridPartType, 2 > >::run(*grid,
      //                                                                                  testcase_1_new_sipdg_2_reference_discretization,
      //                                                                                  num_refinements,
      //                                                                                  testcase_1_exact_solution,
      //                                                                                  "NEW SIP discontinuous
      //                                                                                  galerkin discretization,
      //                                                                                  polOrder 2",
      //                                                                                  "testcase_1");
      //      info << std::endl;
      //      const Example::NewSIPDGDiscretization< GridPartType, 3 >
      //          testcase_1_new_sipdg_3_reference_discretization(*reference_grid_part,
      //                                                          testcase_1_boundary_info,
      //                                                          testcase_1_diffusion,
      //                                                          testcase_1_force,
      //                                                          testcase_1_dirichlet,
      //                                                          testcase_1_neumann);
      //      ConvergenceStudy< Example::NewSIPDGDiscretization< GridPartType, 3 > >::run(*grid,
      //                                                                                  testcase_1_new_sipdg_3_reference_discretization,
      //                                                                                  num_refinements,
      //                                                                                  testcase_1_exact_solution,
      //                                                                                  "NEW SIP discontinuous
      //                                                                                  galerkin discretization,
      //                                                                                  polOrder 3",
      //                                                                                  "testcase_1");
      //      info << std::endl;
      //      // symmetric weighted interior penalty discontinuous galerkin discretization
      //      const SWIPDG_1_DiscretizationType testcase_1_swipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                     testcase_1_boundary_info,
      //                                                                                     testcase_1_diffusion,
      //                                                                                     testcase_1_force,
      //                                                                                     testcase_1_dirichlet,
      //                                                                                     testcase_1_neumann);
      //      ConvergenceStudy< SWIPDG_1_DiscretizationType >::run(*grid,
      //                                                           testcase_1_swipdg_1_reference_discretization,
      //                                                           num_refinements,
      //                                                           testcase_1_exact_solution,
      //                                                           "SWIP discontinuous galerkin discretization, polOrder
      //                                                           1",
      //                                                           "testcase_1");
      //      info << std::endl;

      info << "+==============================================================+" << std::endl;
      info << "|+============================================================+|" << std::endl;
      info << "||  Testcase 2: smooth data, smooth nonhomogeneous dirichlet  ||" << std::endl;
      info << "||              (see page 858 in Epshteyn, Riviere, 2007)     ||" << std::endl;
      info << "|+------------------------------------------------------------+|" << std::endl;
      info << "||  domain = [0, 1] x [0 , 1]                                 ||" << std::endl;
      info << "||  diffusion = 1                                             ||" << std::endl;
      info << "||  force     = 64 pi^2 (cos(8  pi x) + cos(8 pi y))          ||" << std::endl;
      info << "||  dirichlet = cos(8 pi x) + cos(8 pi y)                     ||" << std::endl;
      info << "||  exact solution = cos(8 pi x) + cos(8 pi y)                ||" << std::endl;
      info << "|+============================================================+|" << std::endl;
      info << "+==============================================================+" << std::endl;
      info << std::endl;
      grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(0, 1, 16));
      grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 1); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      const Stuff::GridboundaryAllDirichlet<typename GridPartType::IntersectionType> testcase_2_boundary_info;
      const ConstantFunctionType testcase_2_diffusion(1.0);
      const ExpressionFunctionType testcase_2_force(
          "x", "64.0 * pi * pi * (cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1]))", integration_order);
      const ExpressionFunctionType testcase_2_dirichlet(
          "x", "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])", integration_order);
      const ConstantFunctionType testcase_2_neumann(1.0);
      const ExpressionFunctionType testcase_2_exact_solution(
          "x",
          "cos(8.0 * pi * x[0]) + cos(8.0 * pi * x[1])",
          {"-8.0 * pi * sin(8.0 * pi * x[0])", "-8.0 * pi * sin(8.0 * pi * x[1])"},
          integration_order);

      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_2_cg_1_reference_discretization(*reference_grid_part,
                                                                             testcase_2_boundary_info,
                                                                             testcase_2_diffusion,
                                                                             testcase_2_force,
                                                                             testcase_2_dirichlet,
                                                                             testcase_2_neumann);
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_2_cg_1_reference_discretization,
                                                     num_refinements,
                                                     testcase_2_exact_solution,
                                                     "continuous galerkin discretization, polOrder 1",
                                                     "testcase_2");
      info << std::endl;
      // symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_1_DiscretizationType testcase_2_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                   testcase_2_boundary_info,
                                                                                   testcase_2_diffusion,
                                                                                   testcase_2_force,
                                                                                   testcase_2_dirichlet,
                                                                                   testcase_2_neumann);
      ConvergenceStudy<SIPDG_1_DiscretizationType>::run(*grid,
                                                        testcase_2_sipdg_1_reference_discretization,
                                                        num_refinements,
                                                        testcase_2_exact_solution,
                                                        "SIP discontinuous galerkin discretization, polOrder 1",
                                                        "testcase_2");
      info << std::endl;
      //      // symmetric weighted interior penalty discontinuous galerkin discretization
      //      const SWIPDG_1_DiscretizationType testcase_2_swipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                     testcase_2_boundary_info,
      //                                                                                     testcase_2_diffusion,
      //                                                                                     testcase_2_force,
      //                                                                                     testcase_2_dirichlet,
      //                                                                                     testcase_2_neumann);
      //      ConvergenceStudy< SWIPDG_1_DiscretizationType >::run(*grid,
      //                                                           testcase_2_swipdg_1_reference_discretization,
      //                                                           num_refinements,
      //                                                           testcase_2_exact_solution,
      //                                                           "SWIP discontinuous galerkin discretization, polOrder
      //                                                           1",
      //                                                           "testcase_2");
      //      info << std::endl;

      info << "+============================================================+" << std::endl;
      info << "|+==========================================================+|" << std::endl;
      info << "||  Testcase 3: constant data, mixed dirichlet and neumann  ||" << std::endl;
      info << "|+----------------------------------------------------------+|" << std::endl;
      info << "||  domain = [0, 1] x [0 , 1]                               ||" << std::endl;
      info << "||  diffusion = 1                                           ||" << std::endl;
      info << "||  force     = 1                                           ||" << std::endl;
      info << "||  neumann   = 0.1       on the right side                 ||" << std::endl;
      info << "||  dirichlet = 1/4 * x*y everywhere else                   ||" << std::endl;
      info << "||  reference solution: CG solution on finest grid          ||" << std::endl;
      info << "|+==========================================================+|" << std::endl;
      info << "+============================================================+" << std::endl;
      info << std::endl;
      grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(0, 1, 2));
      grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 1); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      DomainType neumann_normal(0.0);
      neumann_normal[0] = 1.0;
      const Stuff::GridboundaryNormalBased<typename GridPartType::IntersectionType> testcase_3_boundary_info(
          true, {}, {neumann_normal});
      const ConstantFunctionType testcase_3_diffusion(1.0);
      const ConstantFunctionType testcase_3_force(1.0);
      const ExpressionFunctionType testcase_3_dirichlet("x", "0.25 * x[0] * x[1]", integration_order);
      const ConstantFunctionType testcase_3_neumann(0.1);
      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_3_cg_1_reference_discretization(*reference_grid_part,
                                                                             testcase_3_boundary_info,
                                                                             testcase_3_diffusion,
                                                                             testcase_3_force,
                                                                             testcase_3_dirichlet,
                                                                             testcase_3_neumann);
      const auto testcase_3_reference_solution = testcase_3_cg_1_reference_discretization.solve();
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_3_cg_1_reference_discretization,
                                                     num_refinements,
                                                     *testcase_3_reference_solution,
                                                     "continuous galerkin discretization, polOrder 1",
                                                     "testcase_3");
      info << std::endl;
      // adhoc symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_Adhoc_1_DiscretizationType testcase_3_adhoc_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                               testcase_3_boundary_info,
                                                                                               testcase_3_diffusion,
                                                                                               testcase_3_force,
                                                                                               testcase_3_dirichlet,
                                                                                               testcase_3_neumann);
      ConvergenceStudy<SIPDG_Adhoc_1_DiscretizationType>::run(
          *grid,
          testcase_3_adhoc_sipdg_1_reference_discretization,
          num_refinements,
          *testcase_3_reference_solution,
          "Adhoc SIP discontinuous galerkin discretization, polOrder 1",
          "testcase_3");
      info << std::endl;
      // symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_1_DiscretizationType testcase_3_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                   testcase_3_boundary_info,
                                                                                   testcase_3_diffusion,
                                                                                   testcase_3_force,
                                                                                   testcase_3_dirichlet,
                                                                                   testcase_3_neumann);
      ConvergenceStudy<SIPDG_1_DiscretizationType>::run(*grid,
                                                        testcase_3_sipdg_1_reference_discretization,
                                                        num_refinements,
                                                        *testcase_3_reference_solution,
                                                        "SIP discontinuous galerkin discretization, polOrder 1",
                                                        "testcase_3");
      info << std::endl;
      //      // symmetric weighted interior penalty discontinuous galerkin discretization
      //      const SWIPDG_1_DiscretizationType testcase_3_swipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                     testcase_3_boundary_info,
      //                                                                                     testcase_3_diffusion,
      //                                                                                     testcase_3_force,
      //                                                                                     testcase_3_dirichlet,
      //                                                                                     testcase_3_neumann);
      //      ConvergenceStudy< SWIPDG_1_DiscretizationType >::run(*grid,
      //                                                           testcase_3_swipdg_1_reference_discretization,
      //                                                           num_refinements,
      //                                                           *testcase_3_reference_solution,
      //                                                           "SWIP discontinuous galerkin discretization, polOrder
      //                                                           1",
      //                                                           "testcase_3");
      //      info << std::endl;

      info << "+==================================================================================+" << std::endl;
      info << "|+================================================================================+|" << std::endl;
      info << "||  Testcase 4: local thermal block problem                                       ||" << std::endl;
      info << "||              (see http://wwwmath.uni-muenster.de/num/publications/2013/AO13/)  ||" << std::endl;
      info << "|+--------------------------------------------------------------------------------+|" << std::endl;
      info << "||  domain = [0, 1] x [0 , 1]                                                     ||" << std::endl;
      info << "||  diffusion:  see page 3                                                        ||" << std::endl;
      info << "||  force     = 1                                                                 ||" << std::endl;
      info << "||  dirichlet = 0                                                                 ||" << std::endl;
      info << "||  reference solution: CG solution on finest grid                                ||" << std::endl;
      info << "|+================================================================================+|" << std::endl;
      info << "+==================================================================================+" << std::endl;
      info << std::endl;
      grid_provider = std::unique_ptr<GridProviderType>(new GridProviderType(0, 1, 6));
      grid = grid_provider->grid();
      grid->globalRefine(1);
      for (size_t ii = 1; ii <= (num_refinements + 1); ++ii)
        grid->globalRefine(GridType::refineStepsForHalf);
      reference_grid_part = std::unique_ptr<GridPartType>(new GridPartType(*grid, grid->maxLevel()));
      vtk_writer =
          std::unique_ptr<VTKWriterType>(new VTKWriterType(reference_grid_part->gridView(), VTK::nonconforming));

      const Stuff::GridboundaryAllDirichlet<typename GridPartType::IntersectionType> testcase_4_boundary_info;
      const RangeFieldType mu_zero = 0.1;
      const RangeFieldType mu_one  = 1.0;
      const RangeFieldType mu_two = 0.01;
      const CheckerboardFunctionType testcase_4_diffusion(DomainType(0.0),
                                                          DomainType(1.0),
                                                          {6, 6},
                                                          {mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero,
                                                           mu_two,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_one,
                                                           mu_zero,
                                                           mu_zero});
      const ConstantFunctionType testcase_4_force(1.0);
      const ConstantFunctionType testcase_4_dirichlet(0.0);
      const ConstantFunctionType testcase_4_neumann(0.0);
      // continuous galerkin discretization
      const CG_1_DiscretizationType testcase_4_cg_1_reference_discretization(*reference_grid_part,
                                                                             testcase_4_boundary_info,
                                                                             testcase_4_diffusion,
                                                                             testcase_4_force,
                                                                             testcase_4_dirichlet,
                                                                             testcase_4_neumann);
      const auto testcase_4_reference_solution = testcase_4_cg_1_reference_discretization.solve();
      ConvergenceStudy<CG_1_DiscretizationType>::run(*grid,
                                                     testcase_4_cg_1_reference_discretization,
                                                     num_refinements,
                                                     *testcase_4_reference_solution,
                                                     "continuous galerkin discretization, polOrder 1",
                                                     "testcase_4");
      info << std::endl;
      // adhoc symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_Adhoc_1_DiscretizationType testcase_4_adhoc_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                               testcase_4_boundary_info,
                                                                                               testcase_4_diffusion,
                                                                                               testcase_4_force,
                                                                                               testcase_4_dirichlet,
                                                                                               testcase_4_neumann);
      ConvergenceStudy<SIPDG_Adhoc_1_DiscretizationType>::run(
          *grid,
          testcase_4_adhoc_sipdg_1_reference_discretization,
          num_refinements,
          *testcase_4_reference_solution,
          "Adhoc SIP discontinuous galerkin discretization, polOrder 1",
          "testcase_4");
      info << std::endl;
      // symmetric interior penalty discontinuous galerkin discretization
      const SIPDG_1_DiscretizationType testcase_4_sipdg_1_reference_discretization(*reference_grid_part,
                                                                                   testcase_4_boundary_info,
                                                                                   testcase_4_diffusion,
                                                                                   testcase_4_force,
                                                                                   testcase_4_dirichlet,
                                                                                   testcase_4_neumann);
      ConvergenceStudy<SIPDG_1_DiscretizationType>::run(*grid,
                                                        testcase_4_sipdg_1_reference_discretization,
                                                        num_refinements,
                                                        *testcase_4_reference_solution,
                                                        "SIP discontinuous galerkin discretization, polOrder 1",
                                                        "testcase_4");
      info << std::endl;
      //      // symmetric weighted interior penalty discontinuous galerkin discretization
      //      const SWIPDG_1_DiscretizationType testcase_4_swipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                     testcase_4_boundary_info,
      //                                                                                     testcase_4_diffusion,
      //                                                                                     testcase_4_force,
      //                                                                                     testcase_4_dirichlet,
      //                                                                                     testcase_4_neumann);
      //      ConvergenceStudy< SWIPDG_1_DiscretizationType >::run(*grid,
      //                                                           testcase_4_swipdg_1_reference_discretization,
      //                                                           num_refinements,
      //                                                           *testcase_4_reference_solution,
      //                                                           "SWIP discontinuous galerkin discretization, polOrder
      //                                                           1",
      //                                                           "testcase_4");
      //      info << std::endl;

      //      if (spe_10_model_1_datafile_found) {
      //        info << "+=====================================================================+" << std::endl;
      //        info << "|+===================================================================+|" << std::endl;
      //        info << "||  Testcase 5: SPE10                                                ||" << std::endl;
      //        info << "||              (see http://www.spe.org/web/csp/datasets/set01.htm)  ||" << std::endl;
      //        info << "|+-------------------------------------------------------------------+|" << std::endl;
      //        info << "||  domain = [0, 5] x [0 , 1]                                        ||" << std::endl;
      //        info << "||  diffusion:  spe10 model 1                                        ||" << std::endl;
      //        info << "||  force     = 1                                                    ||" << std::endl;
      //        info << "||  dirichlet = 0                                                    ||" << std::endl;
      //        info << "||  reference solution: CG solution on finest grid                   ||" << std::endl;
      //        info << "|+===================================================================+|" << std::endl;
      //        info << "+=====================================================================+" << std::endl;
      //        DomainType upper_right(1.0);
      //        upper_right[0] = 5.0;
      //        grid_provider = std::unique_ptr< GridProviderType >(new GridProviderType(DomainType(0.0),
      //                                                                                 upper_right,
      //                                                                                 {100u, 20u}));
      //        grid = grid_provider->grid();
      //        grid->globalRefine(1);
      //        for (size_t ii = 1; ii <= (num_refinements + 1); ++ii)
      //          grid->globalRefine(GridType::refineStepsForHalf);
      //        reference_grid_part = std::unique_ptr< GridPartType >(new GridPartType(*grid, grid->maxLevel()));
      //        vtk_writer = std::unique_ptr< VTKWriterType >(new VTKWriterType(reference_grid_part->gridView(),
      //                                                                        VTK::nonconforming));

      //        const Stuff::GridboundaryAllDirichlet< typename GridPartType::IntersectionType >
      //        testcase_5_boundary_info;
      //        const SPE10Model1FunctionType   testcase_5_diffusion(spe_10_model_1_datafile,
      //                                                             DomainType(0.0), upper_right);
      //        const ExpressionFunctionType    testcase_5_force("x",
      //                                                         "100.0*exp(-1.0*((((x[0]-0.95)*(x[0]-0.95))+((x[1]-0.65)*(x[1]-0.65)))/(2*0.05*0.05)))-100.0*exp(-1.0*((((x[0]-4.3)*(x[0]-4.3))+((x[1]-0.35)*(x[1]-0.35)))/(2*0.05*0.05)))",
      //                                                         integration_order);
      //        const ConstantFunctionType      testcase_5_dirichlet(0.0);
      //        const ConstantFunctionType      testcase_5_neumann(0.0);
      //        // continuous galerkin discretization
      //        const CG_1_DiscretizationType testcase_5_cg_1_reference_discretization(*reference_grid_part,
      //                                                                               testcase_5_boundary_info,
      //                                                                               testcase_5_diffusion,
      //                                                                               testcase_5_force,
      //                                                                               testcase_5_dirichlet,
      //                                                                               testcase_5_neumann);
      //        const auto testcase_5_reference_solution = testcase_5_cg_1_reference_discretization.solve();
      //        ConvergenceStudy< CG_1_DiscretizationType >::run(*grid,
      //                                                         testcase_5_cg_1_reference_discretization,
      //                                                         num_refinements,
      //                                                         *testcase_5_reference_solution,
      //                                                         "continuous galerkin discretization, polOrder 1",
      //                                                         "testcase_5");
      //        info << std::endl;
      //        // symmetric interior penalty discontinuous galerkin discretization
      //        const SIPDG_1_DiscretizationType testcase_5_sipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                     testcase_5_boundary_info,
      //                                                                                     testcase_5_diffusion,
      //                                                                                     testcase_5_force,
      //                                                                                     testcase_5_dirichlet,
      //                                                                                     testcase_5_neumann);
      //        ConvergenceStudy< SIPDG_1_DiscretizationType >::run(*grid,
      //                                                            testcase_5_sipdg_1_reference_discretization,
      //                                                            num_refinements,
      //                                                            *testcase_5_reference_solution,
      //                                                            "SIP discontinuous galerkin discretization, polOrder
      //                                                            1",
      //                                                            "testcase_5");
      //        info << std::endl;
      //        // symmetric weighted interior penalty discontinuous galerkin discretization
      //        const SWIPDG_1_DiscretizationType testcase_5_swipdg_1_reference_discretization(*reference_grid_part,
      //                                                                                       testcase_5_boundary_info,
      //                                                                                       testcase_5_diffusion,
      //                                                                                       testcase_5_force,
      //                                                                                       testcase_5_dirichlet,
      //                                                                                       testcase_5_neumann);
      //        ConvergenceStudy< SWIPDG_1_DiscretizationType >::run(*grid,
      //                                                             testcase_5_swipdg_1_reference_discretization,
      //                                                             num_refinements,
      //                                                             *testcase_5_reference_solution,
      //                                                             "SWIP discontinuous galerkin discretization,
      //                                                             polOrder 1",
      //                                                             "testcase_5");
      //        info << std::endl;
      //      } // spe 10 testcase
    } // read or write settings file

    // done
    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
    std::abort();
  } catch (std::exception& e) {
    std::cerr << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    std::abort();
  } // try
} // main
