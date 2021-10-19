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

    // output
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);
    std::string filename = config.get("output.prefix", "cpp");
    filename += XT::Grid::is_yaspgrid<typename CellModelSolver::G>::value ? "_cube" : "_simplex";
    filename += XT::Common::to_string(num_elements_x) + "x" + XT::Common::to_string(num_elements_y);
    filename += "_Be" + XT::Common::to_string(Be);
    filename += "_Ca" + XT::Common::to_string(Ca);
    filename += "_Pa" + XT::Common::to_string(Pa);
    filename += "_Fa" + XT::Common::to_string(Fa);
    filename += "_dt" + XT::Common::to_string(dt);
    filename += "_writestep" + XT::Common::to_string(write_step);
    filename += "_polorder" + XT::Common::to_string(pol_order);

    // solver
    const double gmres_reduction = DXTC_CONFIG_GET("gmres_reduction", 1e-13);
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

      // set ofield and stokes values
      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_pfield_dofs(pfield_output_dofs, kk);
        const auto& pfield_source_dofs = model_solver.pfield_deim_source_dofs(kk)[0];
        const auto& pfield_ofield_source_dofs = model_solver.pfield_deim_source_dofs(kk)[1];
        const auto& pfield_stokes_source_dofs = model_solver.pfield_deim_source_dofs(kk)[2];
        std::vector<double> pfield_restricted_source(pfield_source_dofs.size());
        std::vector<double> pfield_ofield_restricted_source(pfield_ofield_source_dofs.size());
        std::vector<double> pfield_stokes_restricted_source(pfield_stokes_source_dofs.size());
        for (size_t ll = 0; ll < pfield_source_dofs.size(); ++ll)
          pfield_restricted_source[ll] = pfield_source[pfield_source_dofs[ll]];
        model_solver.set_pfield_vec_dofs(kk, pfield_restricted_source, pfield_source_dofs);
        model_solver2.set_pfield_vec(kk, pfield_source);
        for (size_t ll = 0; ll < pfield_ofield_source_dofs.size(); ++ll)
          pfield_ofield_restricted_source[ll] = ofield_source[pfield_ofield_source_dofs[ll]];
        model_solver.set_ofield_vec_dofs(kk, pfield_ofield_restricted_source, pfield_ofield_source_dofs);
        model_solver2.set_ofield_vec(kk, ofield_source);
        for (size_t ll = 0; ll < pfield_stokes_source_dofs.size(); ++ll)
          pfield_stokes_restricted_source[ll] = stokes_source[pfield_stokes_source_dofs[ll]];
        model_solver.set_stokes_vec_dofs(pfield_stokes_restricted_source, pfield_stokes_source_dofs);
        model_solver2.set_stokes_vec(stokes_source);
      }

      // change source_vectors randomly to avoid masking errors
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

      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_pfield_dofs(pfield_output_dofs, kk);
        auto begin = std::chrono::steady_clock::now();
        model_solver.prepare_pfield_operator(kk, true);
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
        begin = std::chrono::steady_clock::now();
        model_solver2.prepare_pfield_operator(kk, false);
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
      // set ofield and stokes values
      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_ofield_dofs(ofield_output_dofs, kk);
        const auto& ofield_pfield_source_dofs = model_solver.ofield_deim_source_dofs(kk)[0];
        const auto& ofield_source_dofs = model_solver.ofield_deim_source_dofs(kk)[1];
        const auto& ofield_stokes_source_dofs = model_solver.ofield_deim_source_dofs(kk)[2];
        std::vector<double> ofield_pfield_restricted_source(ofield_pfield_source_dofs.size());
        std::vector<double> ofield_restricted_source(ofield_source_dofs.size());
        std::vector<double> ofield_stokes_restricted_source(ofield_stokes_source_dofs.size());
        for (size_t ll = 0; ll < ofield_pfield_source_dofs.size(); ++ll)
          ofield_pfield_restricted_source[ll] = pfield_source[ofield_pfield_source_dofs[ll]];
        model_solver.set_pfield_vec_dofs(kk, ofield_pfield_restricted_source, ofield_pfield_source_dofs);
        model_solver2.set_pfield_vec(kk, pfield_source);
        for (size_t ll = 0; ll < ofield_source_dofs.size(); ++ll)
          ofield_restricted_source[ll] = ofield_source[ofield_source_dofs[ll]];
        model_solver.set_ofield_vec_dofs(kk, ofield_restricted_source, ofield_source_dofs);
        model_solver2.set_ofield_vec(kk, ofield_source);
        for (size_t ll = 0; ll < ofield_stokes_source_dofs.size(); ++ll)
          ofield_stokes_restricted_source[ll] = stokes_source[ofield_stokes_source_dofs[ll]];
        model_solver.set_stokes_vec_dofs(ofield_stokes_restricted_source, ofield_stokes_source_dofs);
        model_solver2.set_stokes_vec(stokes_source);
      }

      // change source_vectors randomly to avoid masking errors
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

      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.compute_restricted_ofield_dofs(ofield_output_dofs, kk);
        auto begin = std::chrono::steady_clock::now();
        model_solver.prepare_ofield_operator(kk, true);
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
        model_solver2.prepare_ofield_operator(kk, false);
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
      }
      model_solver.compute_restricted_stokes_dofs(stokes_output_dofs);
      const auto& stokes_pfield_source_dofs = model_solver.stokes_deim_source_dofs()[0];
      const auto& stokes_ofield_source_dofs = model_solver.stokes_deim_source_dofs()[1];
      const auto& stokes_source_dofs = model_solver.stokes_deim_source_dofs()[2];
      std::vector<double> stokes_pfield_restricted_source(stokes_pfield_source_dofs.size());
      std::vector<double> stokes_ofield_restricted_source(stokes_ofield_source_dofs.size());
      // set pfield and ofield values
      for (size_t kk = 0; kk < num_cells; ++kk) {
        for (size_t ll = 0; ll < stokes_pfield_source_dofs.size(); ++ll)
          stokes_pfield_restricted_source[ll] = pfield_source[stokes_pfield_source_dofs[ll]];
        model_solver.set_pfield_vec_dofs(kk, stokes_pfield_restricted_source, stokes_pfield_source_dofs);
        model_solver2.set_pfield_vec(kk, pfield_source);
        for (size_t ll = 0; ll < stokes_ofield_source_dofs.size(); ++ll)
          stokes_ofield_restricted_source[ll] = ofield_source[stokes_ofield_source_dofs[ll]];
        model_solver.set_ofield_vec_dofs(kk, stokes_ofield_restricted_source, stokes_ofield_source_dofs);
        model_solver2.set_ofield_vec(kk, ofield_source);
      }

      // change source_vectors randomly to avoid masking errors
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

      std::this_thread::sleep_for(std::chrono::seconds(2));
      auto begin = std::chrono::steady_clock::now();
      VectorType restricted_result, restricted_jac_result;
      for (size_t ii = 0; ii < 10000; ii++) {
        model_solver.prepare_stokes_operator(true);
        stokes_restricted_prep_time += std::chrono::steady_clock::now() - begin;
        const size_t num_source_dofs = stokes_source_dofs.size();
        VectorType restricted_source(num_source_dofs, 0.);
        for (size_t jj = 0; jj < num_source_dofs; ++jj)
          restricted_source[jj] = stokes_source[stokes_source_dofs[jj]];
        begin = std::chrono::steady_clock::now();
        restricted_result = model_solver.apply_stokes_operator(restricted_source, true);
        stokes_restricted_apply_time += std::chrono::steady_clock::now() - begin;
        begin = std::chrono::steady_clock::now();
        restricted_jac_result = model_solver.apply_stokes_jacobian(restricted_source, true);
        stokes_restricted_jac_time += std::chrono::steady_clock::now() - begin;
      }
      std::this_thread::sleep_for(std::chrono::seconds(2));
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
    std::cout << "Num dofs: "
              << "Pfield: " << pfield_size << ", Ofield: " << ofield_size << ", Stokes: " << stokes_size << std::endl;
    std::cout << "Pfield:" << std::endl;
    const auto& pfield_source_dofs = model_solver.pfield_deim_source_dofs(0);
    const size_t pfield_pfield_constrained_size = pfield_source_dofs[0].size();
    const size_t pfield_ofield_constrained_size = pfield_source_dofs[1].size();
    const size_t pfield_stokes_constrained_size = pfield_source_dofs[2].size();
    std::cout << "Constrained input dofs: "
              << "Pfield: " << pfield_pfield_constrained_size << " (factor "
              << double(pfield_size) / double(pfield_pfield_constrained_size) << "), "
              << "Ofield: " << pfield_ofield_constrained_size << " (factor "
              << double(ofield_size) / double(pfield_ofield_constrained_size) << "), "
              << "Stokes: " << pfield_stokes_constrained_size << " (factor "
              << double(stokes_size) / double(pfield_stokes_constrained_size) << "), " << std::endl;
    std::cout << "prep: " << pfield_prep_time.count() << "  vs. " << pfield_restricted_prep_time.count()
              << ", factor: " << pfield_prep_time.count() / pfield_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << pfield_apply_time.count() << "  vs. " << pfield_restricted_apply_time.count()
              << ", factor: " << pfield_apply_time.count() / pfield_restricted_apply_time.count() << std::endl;
    std::cout << "jac: " << pfield_jac_time.count() << "  vs. " << pfield_restricted_jac_time.count()
              << ", factor: " << pfield_jac_time.count() / pfield_restricted_jac_time.count() << std::endl;
    std::cout << "Ofield:" << std::endl;
    const auto& ofield_source_dofs = model_solver.ofield_deim_source_dofs(0);
    const size_t ofield_pfield_constrained_size = ofield_source_dofs[0].size();
    const size_t ofield_ofield_constrained_size = ofield_source_dofs[1].size();
    const size_t ofield_stokes_constrained_size = ofield_source_dofs[2].size();
    std::cout << "Constrained input dofs: "
              << "Pfield: " << ofield_pfield_constrained_size << " (factor "
              << double(pfield_size) / double(ofield_pfield_constrained_size) << "), "
              << "Ofield: " << ofield_ofield_constrained_size << " (factor "
              << double(ofield_size) / double(ofield_ofield_constrained_size) << "), "
              << "Stokes: " << ofield_stokes_constrained_size << " (factor "
              << double(stokes_size) / double(ofield_stokes_constrained_size) << "), " << std::endl;
    std::cout << "prep: " << ofield_prep_time.count() << "  vs. " << ofield_restricted_prep_time.count()
              << ", factor: " << ofield_prep_time.count() / ofield_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << ofield_apply_time.count() << "  vs. " << ofield_restricted_apply_time.count()
              << ", factor: " << ofield_apply_time.count() / ofield_restricted_apply_time.count() << std::endl;
    std::cout << "jac: " << ofield_jac_time.count() << "  vs. " << ofield_restricted_jac_time.count()
              << ", factor: " << ofield_jac_time.count() / ofield_restricted_jac_time.count() << std::endl;
    std::cout << "Stokes:" << std::endl;
    const auto& stokes_source_dofs = model_solver.stokes_deim_source_dofs();
    const size_t stokes_pfield_constrained_size = stokes_source_dofs[0].size();
    const size_t stokes_ofield_constrained_size = stokes_source_dofs[1].size();
    const size_t stokes_stokes_constrained_size = stokes_source_dofs[2].size();
    std::cout << "Constrained input dofs: "
              << "Pfield: " << stokes_pfield_constrained_size << " (factor "
              << double(pfield_size) / double(stokes_pfield_constrained_size) << "), "
              << "Ofield: " << stokes_ofield_constrained_size << " (factor "
              << double(ofield_size) / double(stokes_ofield_constrained_size) << "), "
              << "Stokes: " << stokes_stokes_constrained_size << " (factor "
              << double(stokes_size) / double(stokes_stokes_constrained_size) << "), " << std::endl;
    std::cout << "prep: " << stokes_prep_time.count() << "  vs. " << stokes_restricted_prep_time.count()
              << ", factor: " << stokes_prep_time.count() / stokes_restricted_prep_time.count() << std::endl;
    std::cout << "apply: " << stokes_apply_time.count() << "  vs. " << stokes_restricted_apply_time.count()
              << ", factor: " << stokes_apply_time.count() / stokes_restricted_apply_time.count() << std::endl;
    std::cout << "jac: " << stokes_jac_time.count() << "  vs. " << stokes_restricted_jac_time.count()
              << ", factor: " << stokes_jac_time.count() / stokes_restricted_jac_time.count() << std::endl;
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
