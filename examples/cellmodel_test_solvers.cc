// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <dune/gdt/test/cellmodel/cellmodel.hh>

#include <chrono>
#include <random>

int main(int argc, char* argv[])
{
  using namespace Dune;
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    const size_t num_threads = DXTC_CONFIG.get("threading.max_count", 1u);
#if HAVE_TBB
    DXTC_CONFIG.set("threading.partition_factor", 1, true);
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
    double t_end = config.template get<double>("fem.t_end", 0.1);
    double dt = config.template get<double>("fem.dt", 0.01);
    const int pol_order = config.template get<int>("fem.degree", 1, 0, 0);
    const int overintegrate = config.template get<int>("fem.superintegration_order", 2, 0, 0);

    // problem parameters
    double L = config.template get<double>("problem.L", 1e-6);
    double U = config.template get<double>("problem.U", 1e-6);
    double rho = config.template get<double>("problem.rho", 1.e3);
    double eta = config.template get<double>("problem.eta", 2.e3);
    double sigma = config.template get<double>("problem.sigma", 0.0188);
    double b_N = config.template get<double>("problem.b_N", 1.26e-14);
    double k = config.template get<double>("problem.k", 2.e-9);
    double xi = config.template get<double>("problem.xi", 1.1);
    double eta_rot = config.template get<double>("problem.eta_rot", 3.3e3);
    double zeta = config.template get<double>("problem.zeta", 2.e3);
    double epsilon = config.template get<double>("problem.epsilon", 0.21);
    double gamma = config.template get<double>("problem.gamma", 0.025);
    double c_1 = config.template get<double>("problem.c_1", 5.);
    double Re = rho * U * L / eta;
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

    // output
    std::string filename = config.get("output.filename", "cellmodel");
    bool subsampling = config.get<bool>("output.subsampling", true);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);

    // GMRES options
    const double gmres_reduction = DXTC_CONFIG_GET("gmres_reduction", 1e-10);
    const int gmres_restart = DXTC_CONFIG_GET("gmres_restart", 100);
    const double inner_gmres_reduction = DXTC_CONFIG_GET("inner_gmres_reduction", 1e-3);
    const int inner_gmres_maxit = DXTC_CONFIG_GET("inner_gmres_maxit", 10);
    const int gmres_verbose = DXTC_CONFIG_GET("gmres_verbose", 0);
    const CellModelLinearSolverType pfield_solver_type =
        string_to_solver_type(DXTC_CONFIG_GET("pfield_solver_type", "gmres"));
    const CellModelMassMatrixSolverType pfield_mass_matrix_solver_type =
        string_to_mass_matrix_solver_type(DXTC_CONFIG_GET("pfield_mass_matrix_solver_type", "sparse_ldlt"));
    const CellModelLinearSolverType ofield_solver_type =
        string_to_solver_type(DXTC_CONFIG_GET("ofield_solver_type", "schur_gmres"));
    const CellModelMassMatrixSolverType ofield_mass_matrix_solver_type =
        string_to_mass_matrix_solver_type(DXTC_CONFIG_GET("ofield_mass_matrix_solver_type", "sparse_ldlt"));
    const std::string stokes_solver_type_string = DXTC_CONFIG_GET("stokes_solver_type", "schur_cg_A_direct_prec_mass");
    const StokesSolverType stokes_solver_type = string_to_stokes_solver_type(stokes_solver_type_string);

    CellModelSolver model_solver(testcase,
                                 t_end,
                                 dt,
                                 num_elements_x,
                                 num_elements_y,
                                 pol_order,
                                 false,
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
                                 CellModelLinearSolverType::direct,
                                 CellModelMassMatrixSolverType::sparse_ldlt,
                                 CellModelLinearSolverType::direct,
                                 CellModelMassMatrixSolverType::sparse_ldlt,
                                 StokesSolverType::direct,
                                 gmres_reduction,
                                 gmres_restart,
                                 gmres_verbose,
                                 inner_gmres_reduction,
                                 inner_gmres_maxit,
                                 gmres_verbose);
    CellModelSolver model_solver2(testcase,
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
    std::chrono::duration<double> ref_time(0.), time(0.), time_after_reset(0.);
    auto begin = std::chrono::steady_clock::now();
    const auto result2 = model_solver2.solve(false, write_step, filename, subsampling);
    time = std::chrono::steady_clock::now() - begin;
    model_solver2.reset();
    begin = std::chrono::steady_clock::now();
    const auto result_reset = model_solver2.solve(false, write_step, filename, subsampling);
    time_after_reset = std::chrono::steady_clock::now() - begin;
    begin = std::chrono::steady_clock::now();
    const auto result1 = model_solver.solve(false, write_step, filename, subsampling);
    ref_time = std::chrono::steady_clock::now() - begin;
    std::cout << "Timings: ref =  " << ref_time.count() << "  vs. tested solvers = " << time.count()
              << " vs. tested solvers after reset " << time_after_reset.count() << std::endl;
    for (size_t ii = 0; ii < result1.size(); ++ii) {
      for (size_t jj = 0; jj < result1[ii].size(); ++jj) {
        for (size_t kk = 0; kk < result1[ii][jj].size(); ++kk) {
          if (XT::Common::FloatCmp::ne(result1[ii][jj][kk], result2[ii][jj][kk], 1e-8, 1e-8))
            std::cout << "Entries (" << ii << ", " << jj << ", " << kk << ") with different solvers differ by "
                      << std::abs(result1[ii][jj][kk] - result2[ii][jj][kk]) << ", " << result1[ii][jj][kk] << " vs. "
                      << result1[ii][jj][kk] << std::endl;
          if (XT::Common::FloatCmp::ne(result2[ii][jj][kk], result_reset[ii][jj][kk], 1e-8, 1e-8))
            std::cout << "Entries (" << ii << ", " << jj << ", " << kk << ") after reset differ by "
                      << std::abs(result2[ii][jj][kk] - result_reset[ii][jj][kk]) << ", " << result2[ii][jj][kk]
                      << " vs. " << result_reset[ii][jj][kk] << std::endl;
        } // kk
      } // jj
    } // ii
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
