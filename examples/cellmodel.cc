// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/gdt/test/cellmodel/cellmodel.hh>

int main(int argc, char* argv[])
{
  using namespace Dune;
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    const size_t num_threads = DXTC_CONFIG.get("threading.max_count", 1u);
#if HAVE_TBB
    DXTC_CONFIG.set("threading.partition_factor", 10, true);
    XT::Common::threadManager().set_max_threads(num_threads);
#endif

    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // read configuration
    XT::Common::Configuration config("activepolargels.ini");

    // get testcase
    auto testcase = config.template get<std::string>("problem.testcase");

    // grid config
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(90));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(90));

    // timestepping
    double t_end = config.template get<double>("fem.t_end", 340.);
    double dt = config.template get<double>("fem.dt", 0.005);

    // orders
    const int pol_order = config.template get<int>("fem.degree", 1, 0, 0);
    const int overintegrate = config.template get<int>("fem.superintegration_order", 2, 0, 0);

    // problem parameters
    const double epsilon = config.template get<double>("problem.epsilon", 0.21);
    const double gamma = config.template get<double>("problem.gamma", 0.025);
    const double c_1 = config.template get<double>("problem.c_1", 5.);
    const double kappa = config.template get<double>("problem.kappa", 1.);
    const double xi = config.template get<double>("problem.xi", 1.1);
    const double Re = 1e-13;
    const double Be = config.template get<double>("problem.Be", 1.0);
    const double Ca = config.template get<double>("problem.Ca", 1.0);
    const double Pa = config.template get<double>("problem.Pa", 1.0);
    const double Fa = config.template get<double>("problem.Fa", 1.0);

    // solver
    const double gmres_reduction = DXTC_CONFIG_GET("gmres_reduction", 1e-13);
    const int gmres_restart = DXTC_CONFIG_GET("gmres_restart", 100);
    const double inner_gmres_reduction = DXTC_CONFIG_GET("inner_gmres_reduction", 1e-3);
    const int inner_gmres_maxit = DXTC_CONFIG_GET("inner_gmres_maxit", 10);
    const int gmres_verbose = DXTC_CONFIG_GET("gmres_verbose", 0);
    const std::string pfield_solver_type_string = DXTC_CONFIG_GET("pfield_solver_type", "gmres");
    const CellModelLinearSolverType pfield_solver_type = string_to_solver_type(pfield_solver_type_string);
    const CellModelMassMatrixSolverType pfield_mass_matrix_solver_type =
        string_to_mass_matrix_solver_type(DXTC_CONFIG_GET("pfield_mass_matrix_solver_type", "sparse_ldlt"));
    const std::string ofield_solver_type_string = DXTC_CONFIG_GET("ofield_solver_type", "schur_gmres");
    const CellModelLinearSolverType ofield_solver_type = string_to_solver_type(ofield_solver_type_string);
    const CellModelMassMatrixSolverType ofield_mass_matrix_solver_type =
        string_to_mass_matrix_solver_type(DXTC_CONFIG_GET("ofield_mass_matrix_solver_type", "sparse_ldlt"));
    const std::string stokes_solver_type_string = DXTC_CONFIG_GET("stokes_solver_type", "schur_cg_A_direct_prec_mass");
    const StokesSolverType stokes_solver_type = string_to_stokes_solver_type(stokes_solver_type_string);

    // output
    bool subsampling = config.get<bool>("output.subsampling", false);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);
    std::string filename = config.get("output.prefix", "cpp");
    filename += "_" + testcase;
    filename += "_pfield_" + pfield_solver_type_string;
    filename += "_ofield_" + ofield_solver_type_string;
    filename += "_stokes_" + stokes_solver_type_string;
    filename += XT::Grid::is_yaspgrid<typename CellModelSolver::G>::value ? "_cube" : "_simplex";
    filename += XT::Common::to_string(num_elements_x) + "x" + XT::Common::to_string(num_elements_y);
    filename += "_Be" + XT::Common::to_string(Be);
    filename += "_Ca" + XT::Common::to_string(Ca);
    filename += "_Pa" + XT::Common::to_string(Pa);
    filename += "_Fa" + XT::Common::to_string(Fa);
    filename += "_dt" + XT::Common::to_string(dt);
    filename += "_writestep" + XT::Common::to_string(write_step);
    filename += "_polorder" + XT::Common::to_string(pol_order);


    CellModelSolver model_solver(testcase,
                                 t_end,
                                 dt,
                                 num_elements_x,
                                 num_elements_y,
                                 pol_order,
                                 num_threads > 1,
                                 Be,
                                 Ca,
                                 Pa,
                                 Re,
                                 Fa,
                                 xi,
                                 kappa,
                                 c_1,
                                 gamma,
                                 epsilon,
                                 overintegrate,
                                 pfield_solver_type,
                                 pfield_mass_matrix_solver_type,
                                 ofield_solver_type,
                                 ofield_mass_matrix_solver_type,
                                 stokes_solver_type,
                                 gmres_reduction,
                                 gmres_restart,
                                 gmres_verbose,
                                 inner_gmres_reduction,
                                 inner_gmres_maxit,
                                 gmres_verbose);

    auto begin = std::chrono::steady_clock::now();
    model_solver.solve_without_storing(true, write_step, filename, subsampling);
    const std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving took: " << time.count() << " s." << std::endl;
  } catch (Exception& e) {
    std::cerr << "\nDUNE reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (std::exception& e) {
    std::cerr << "\nstl reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Unknown error occured!" << std::endl;
    return EXIT_FAILURE;
  } // try
  return EXIT_SUCCESS;
} // ... main(...)
