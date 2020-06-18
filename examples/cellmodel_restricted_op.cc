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
#include <thread>

int main(int argc, char* argv[])
{
  using namespace Dune;
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
#if HAVE_TBB
    DXTC_CONFIG.set("threading.partition_factor", 1, true);
    XT::Common::threadManager().set_max_threads(1);
#endif

    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // read configuration
    XT::Common::Configuration config("activepolargels.ini");

    // get testcase
    auto testcase = config.template get<std::string>("problem.testcase");

    // grid config
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(16));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(4));

    // timestepping
    double t_end = config.template get<double>("fem.t_end", 340.);
    double dt = config.template get<double>("fem.dt", 0.005);

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
    double beta = config.template get<double>("problem.beta", 0.);
    double In = config.template get<double>("problem.In", 1.);
    double Re = rho * U * L / eta;
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

    // output
    std::string filename = config.get("output.filename", "cellmodel");

    const double gmres_reduction = 1e-10;
    const int gmres_restart = 50;
    const double inner_gmres_reduction = 1e-3;
    const int inner_gmres_maxit = 10;
    const int gmres_verbose = 0;
    const CellModelLinearSolverType solver_type = CellModelLinearSolverType::schur_fgmres_gmres;
    const CellModelMassMatrixSolverType mass_matrix_solver_type = CellModelMassMatrixSolverType::sparse_ldlt;

    CellModelSolver model_solver(testcase,
                                 t_end,
                                 num_elements_x,
                                 num_elements_y,
                                 2,
                                 false,
                                 Be,
                                 Ca,
                                 Pa,
                                 Re,
                                 Fa,
                                 xi,
                                 kappa,
                                 c_1,
                                 beta,
                                 gamma,
                                 epsilon,
                                 In,
                                 solver_type,
                                 mass_matrix_solver_type,
                                 solver_type,
                                 mass_matrix_solver_type,
                                 gmres_reduction,
                                 gmres_restart,
                                 gmres_verbose,
                                 inner_gmres_reduction,
                                 inner_gmres_maxit,
                                 gmres_verbose);
    CellModelSolver model_solver2(testcase,
                                  t_end,
                                  num_elements_x,
                                  num_elements_y,
                                  2,
                                  false,
                                  Be,
                                  Ca,
                                  Pa,
                                  Re,
                                  Fa,
                                  xi,
                                  kappa,
                                  c_1,
                                  beta,
                                  gamma,
                                  epsilon,
                                  In,
                                  solver_type,
                                  mass_matrix_solver_type,
                                  solver_type,
                                  mass_matrix_solver_type,
                                  gmres_reduction,
                                  gmres_restart,
                                  gmres_verbose,
                                  inner_gmres_reduction,
                                  inner_gmres_maxit,
                                  gmres_verbose);

    const size_t pfield_size = model_solver.pfield_vec(0).size();
    const size_t ofield_size = model_solver.ofield_vec(0).size();
    const size_t stokes_size = model_solver.stokes_vec().size();
    // choose all dofs
    // const size_t num_output_dofs = model_solver.pfield_vec(0).size();
    // std::vector<size_t> pfield_output_dofs(max_output_dofs);
    // for (size_t ii = 0; ii < max_output_dofs; ++ii)
    // pfield_output_dofs[ii] = ii;
    // choose 50 random dofs
    const size_t num_output_dofs = 20;
    std::vector<size_t> pfield_output_dofs(num_output_dofs);
    std::vector<size_t> ofield_output_dofs(num_output_dofs);
    std::vector<size_t> stokes_output_dofs(num_output_dofs);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<typename std::mt19937::result_type> pfield_distrib(0, pfield_size - 1);
    std::uniform_int_distribution<typename std::mt19937::result_type> ofield_distrib(0, ofield_size - 1);
    std::uniform_int_distribution<typename std::mt19937::result_type> stokes_distrib(0, stokes_size - 1);
    std::uniform_real_distribution<double> double_distrib(-1., 1.);
    using VectorType = typename CellModelSolver::VectorType;
    const size_t num_cells = model_solver.num_cells();
    auto pfield_source = model_solver.pfield_vec(0);
    auto ofield_source = model_solver.ofield_vec(0);
    auto stokes_source = model_solver.stokes_vec();
    auto pfield_state = model_solver.pfield_vec(0);
    auto ofield_state = model_solver.ofield_vec(0);
    // ensure source_vectors are non-zero to avoid masking errors
    for (auto& entry : pfield_source)
      entry += double_distrib(rng);
    for (auto& entry : ofield_source)
      entry += double_distrib(rng);
    for (auto& entry : stokes_source)
      entry += double_distrib(rng);
    for (auto& entry : pfield_state)
      entry += double_distrib(rng);
    for (auto& entry : ofield_state)
      entry += double_distrib(rng);
    std::chrono::duration<double> pfield_restricted_prep_time(0.);
    std::chrono::duration<double> pfield_restricted_apply_time(0.);
    std::chrono::duration<double> pfield_restricted_jac_time(0.);
    std::chrono::duration<double> pfield_prep_time(0.);
    std::chrono::duration<double> pfield_apply_time(0.);
    std::chrono::duration<double> pfield_jac_time(0.);
    std::chrono::duration<double> ofield_restricted_prep_time(0.);
    std::chrono::duration<double> ofield_restricted_apply_time(0.);
    std::chrono::duration<double> ofield_restricted_jac_time(0.);
    std::chrono::duration<double> ofield_prep_time(0.);
    std::chrono::duration<double> ofield_apply_time(0.);
    std::chrono::duration<double> ofield_jac_time(0.);
    std::chrono::duration<double> stokes_restricted_prep_time(0.);
    std::chrono::duration<double> stokes_restricted_apply_time(0.);
    std::chrono::duration<double> stokes_restricted_jac_time(0.);
    std::chrono::duration<double> stokes_prep_time(0.);
    std::chrono::duration<double> stokes_apply_time(0.);
    std::chrono::duration<double> stokes_jac_time(0.);
    for (size_t run = 0; run < 10; ++run) {
      std::cout << "Pfield run " << run << std::endl;
      // std::cout << "Output_dofs: (";
      for (size_t ii = 0; ii < num_output_dofs; ++ii) {
        pfield_output_dofs[ii] = pfield_distrib(rng);
        // std::cout << pfield_output_dofs[ii] << (ii == num_output_dofs - 1 ? ")\n" : ", ");
      }
      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_pfield_dofs(pfield_output_dofs, kk);
        std::cout << "Restricted start" << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(2));
        auto begin = std::chrono::steady_clock::now();
        model_solver.prepare_pfield_operator(dt, kk, true);
        model_solver.set_pfield_jacobian_state(pfield_state, kk, true);
        pfield_restricted_prep_time += std::chrono::steady_clock::now() - begin;
        const auto& pfield_source_dofs = model_solver.pfield_deim_source_dofs(kk)[0];
        const size_t num_source_dofs = pfield_source_dofs.size();
        VectorType restricted_source(num_source_dofs, 0.);
        for (size_t ii = 0; ii < num_source_dofs; ++ii)
          restricted_source[ii] = pfield_source[pfield_source_dofs[ii]];
        begin = std::chrono::steady_clock::now();
        auto restricted_result = model_solver.apply_pfield_operator(restricted_source, kk, true);
        pfield_restricted_apply_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto restricted_jac_result = model_solver.apply_pfield_jacobian(restricted_source, kk, true);
        pfield_restricted_jac_time += std::chrono::steady_clock::now() - begin;
        std::cout << "Restricted end" << std::endl;
        std::this_thread::sleep_for(std::chrono::seconds(2));
        begin = std::chrono::steady_clock::now();
        model_solver2.prepare_pfield_operator(dt, kk, false);
        model_solver2.set_pfield_jacobian_state(pfield_state, kk, false);
        pfield_prep_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto result = model_solver2.apply_pfield_operator(pfield_source, kk, false);
        pfield_apply_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto jac_result = model_solver2.apply_pfield_jacobian(pfield_source, kk, false);
        pfield_jac_time += std::chrono::steady_clock::now() - begin;
        // There are differences of about 1e-13, caused by the different mv methods in assemble_pfield_rhs (mv from
        // Eigen vs mv_restricted);
        const double apply_tol = 1e-12;
        // For apply_pfield_jacobian, the differences are smaller
        const double jac_tol = 1e-14;
        for (size_t ii = 0; ii < num_output_dofs; ++ii) {
          if (XT::Common::FloatCmp::ne(restricted_result[ii], result[pfield_output_dofs[ii]], apply_tol, apply_tol))
            std::cout << "Failed apply restricted: " << ii << ", " << pfield_output_dofs[ii] << ", "
                      << result[pfield_output_dofs[ii]] << ", " << restricted_result[ii] << std::endl;
          if (XT::Common::FloatCmp::ne(restricted_jac_result[ii], jac_result[pfield_output_dofs[ii]], jac_tol, jac_tol))
            std::cout << "Failed apply restricted jacobian: " << ii << ", " << pfield_output_dofs[ii] << ", "
                      << jac_result[pfield_output_dofs[ii]] << ", " << restricted_jac_result[ii] << std::endl;
        } // ii
      } // kk
      std::cout << "Ofield run " << run << std::endl;
      // std::cout << "Ofield output_dofs: (";
      for (size_t ii = 0; ii < num_output_dofs; ++ii) {
        ofield_output_dofs[ii] = ofield_distrib(rng);
        // std::cout << ofield_output_dofs[ii] << (ii == num_output_dofs - 1 ? ")\n" : ", ");
      }
      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_ofield_dofs(ofield_output_dofs, kk);
        auto begin = std::chrono::steady_clock::now();
        model_solver.prepare_ofield_operator(dt, kk, true);
        model_solver.set_ofield_jacobian_state(ofield_state, kk, true);
        ofield_restricted_prep_time += std::chrono::steady_clock::now() - begin;
        const auto& ofield_source_dofs = model_solver.ofield_deim_source_dofs(kk)[1];
        const size_t num_source_dofs = ofield_source_dofs.size();
        VectorType restricted_source(num_source_dofs, 0.);
        for (size_t ii = 0; ii < num_source_dofs; ++ii)
          restricted_source[ii] = ofield_source[ofield_source_dofs[ii]];
        begin = std::chrono::steady_clock::now();
        auto restricted_result = model_solver.apply_ofield_operator(restricted_source, kk, true);
        ofield_restricted_apply_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto restricted_jac_result = model_solver.apply_ofield_jacobian(restricted_source, kk, true);
        ofield_restricted_jac_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        model_solver2.prepare_ofield_operator(dt, kk, false);
        model_solver2.set_ofield_jacobian_state(ofield_state, kk, false);
        ofield_prep_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto result = model_solver2.apply_ofield_operator(ofield_source, kk, false);
        ofield_apply_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        auto jac_result = model_solver2.apply_ofield_jacobian(ofield_source, kk, false);
        ofield_jac_time += std::chrono::steady_clock::now() - begin;
        // There are differences of about 1e-13, caused by the different mv methods in assemble_ofield_rhs (mv from
        // Eigen vs mv_restricted);
        const double apply_tol = 1e-12;
        // For apply_ofield_jacobian, the differences are smaller
        const double jac_tol = 1e-14;
        for (size_t ii = 0; ii < num_output_dofs; ++ii) {
          if (XT::Common::FloatCmp::ne(restricted_result[ii], result[ofield_output_dofs[ii]], apply_tol, apply_tol))
            std::cout << "Failed apply restricted: " << ii << ", " << ofield_output_dofs[ii] << ", "
                      << result[ofield_output_dofs[ii]] << ", " << restricted_result[ii] << std::endl;
          if (XT::Common::FloatCmp::ne(restricted_jac_result[ii], jac_result[ofield_output_dofs[ii]], jac_tol, jac_tol))
            std::cout << "Failed apply restricted jacobian: " << ii << ", " << ofield_output_dofs[ii] << ", "
                      << jac_result[ofield_output_dofs[ii]] << ", " << restricted_jac_result[ii] << std::endl;
        } // ii
      } // kk
      std::cout << "Stokes run " << run << std::endl;
      // std::cout << "Stokes output_dofs: (";
      for (size_t ii = 0; ii < num_output_dofs; ++ii) {
        stokes_output_dofs[ii] = stokes_distrib(rng);
        // std::cout << ofield_output_dofs[ii] << (ii == num_output_dofs - 1 ? ")\n" : ", ");
      }
      model_solver.compute_restricted_stokes_dofs(stokes_output_dofs);
      auto begin = std::chrono::steady_clock::now();
      model_solver.prepare_stokes_operator(true);
      stokes_restricted_prep_time += std::chrono::steady_clock::now() - begin;
      const auto& stokes_source_dofs = model_solver.stokes_deim_source_dofs()[2];
      const size_t num_source_dofs = stokes_source_dofs.size();
      VectorType restricted_source(num_source_dofs, 0.);
      for (size_t ii = 0; ii < num_source_dofs; ++ii)
        restricted_source[ii] = stokes_source[stokes_source_dofs[ii]];
      begin = std::chrono::steady_clock::now();
      auto restricted_result = model_solver.apply_stokes_operator(restricted_source, true);
      stokes_restricted_apply_time += std::chrono::steady_clock::now() - begin;
      begin = std::chrono::steady_clock::now();
      auto restricted_jac_result = model_solver.apply_stokes_jacobian(restricted_source, true);
      stokes_restricted_jac_time += std::chrono::steady_clock::now() - begin;
      begin = std::chrono::steady_clock::now();
      model_solver2.prepare_stokes_operator(false);
      stokes_prep_time += std::chrono::steady_clock::now() - begin;
      begin = std::chrono::steady_clock::now();
      auto result = model_solver2.apply_stokes_operator(stokes_source, false);
      stokes_apply_time += std::chrono::steady_clock::now() - begin;
      begin = std::chrono::steady_clock::now();
      auto jac_result = model_solver2.apply_stokes_jacobian(stokes_source, false);
      stokes_jac_time += std::chrono::steady_clock::now() - begin;
      // There are differences of about 1e-13, caused by the different mv methods in assemble_stokes_rhs (mv from
      // Eigen vs mv_restricted);
      const double apply_tol = 1e-12;
      // For apply_stokes_jacobian, the differences are smaller
      const double jac_tol = 1e-14;
      for (size_t ii = 0; ii < num_output_dofs; ++ii) {
        if (XT::Common::FloatCmp::ne(restricted_result[ii], result[stokes_output_dofs[ii]], apply_tol, apply_tol))
          std::cout << "Failed apply restricted: " << ii << ", " << stokes_output_dofs[ii] << ", "
                    << result[stokes_output_dofs[ii]] << ", " << restricted_result[ii] << std::endl;
        if (XT::Common::FloatCmp::ne(restricted_jac_result[ii], jac_result[stokes_output_dofs[ii]], jac_tol, jac_tol))
          std::cout << "Failed apply restricted jacobian: " << ii << ", " << stokes_output_dofs[ii] << ", "
                    << jac_result[stokes_output_dofs[ii]] << ", " << restricted_jac_result[ii] << std::endl;
      } // ii
    } // runs
    std::cout << "Pfield:" << std::endl;
    std::cout << "prep: " << pfield_prep_time.count() << "  vs. " << pfield_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << pfield_apply_time.count() << "  vs. " << pfield_restricted_apply_time.count()
              << std::endl;
    std::cout << "jac: " << pfield_jac_time.count() << "  vs. " << pfield_restricted_jac_time.count() << std::endl;
    std::cout << "Ofield:" << std::endl;
    std::cout << "prep: " << ofield_prep_time.count() << "  vs. " << ofield_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << ofield_apply_time.count() << "  vs. " << ofield_restricted_apply_time.count()
              << std::endl;
    std::cout << "jac: " << ofield_jac_time.count() << "  vs. " << ofield_restricted_jac_time.count() << std::endl;
    std::cout << "Stokes:" << std::endl;
    std::cout << "prep: " << stokes_prep_time.count() << "  vs. " << stokes_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << stokes_apply_time.count() << "  vs. " << stokes_restricted_apply_time.count()
              << std::endl;
    std::cout << "jac: " << stokes_jac_time.count() << "  vs. " << stokes_restricted_jac_time.count() << std::endl;
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
