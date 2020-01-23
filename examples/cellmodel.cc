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
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(16));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(4));

    // timestepping
    double t_end = config.template get<double>("fem.t_end", 340.);
    double dt = config.template get<double>("fem.dt", 0.005);
    const bool linearize = config.template get<bool>("problem.linearize", false);
    std::cout << "linearize: " << linearize << std::endl;

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
    // double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    // double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    // double Pa = eta * U * L / k;
    double Be = config.template get<double>("problem.Be", 0.3);
    double Ca = config.template get<double>("problem.Ca", 0.1);
    double Pa = config.template get<double>("problem.Pa", 1.0);
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

    // output
    std::string filename = config.get("output.filename", "cellmodel") + (linearize ? "_linearized" : "");
    bool subsampling = config.get<bool>("output.subsampling", true);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);
    const double gmres_reduction = DXTC_CONFIG_GET("gmres_reduction", 1e-10);
    const int gmres_restart = DXTC_CONFIG_GET("gmres_restart", 100);
    const double inner_gmres_reduction = DXTC_CONFIG_GET("inner_gmres_reduction", 1e-3);
    const int inner_gmres_maxit = DXTC_CONFIG_GET("inner_gmres_maxit", 10);
    const int gmres_verbose = DXTC_CONFIG_GET("gmres_verbose", 0);
    // const std::string pfield_solver_type = DXTC_CONFIG_GET("pfield_solver_type", "custom");

    CellModelSolver model_solver(testcase,
                                 t_end,
                                 num_elements_x,
                                 num_elements_y,
                                 num_threads > 1,
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
                                 PfieldLinearSolverType::schur_fgmres_gmres,
                                 PfieldMassMatrixSolverType::sparse_lu,
                                 "schur",
                                 gmres_reduction,
                                 gmres_restart,
                                 gmres_verbose,
                                 inner_gmres_reduction,
                                 inner_gmres_maxit,
                                 gmres_verbose,
                                 linearize);
#if 1
    auto begin = std::chrono::steady_clock::now();
    auto result = model_solver.solve(dt, true, write_step, filename, subsampling);
    const std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving took: " << time.count() << " s." << std::endl;
    for (const auto& vec1 : result[0])
      std::cout << &vec1 << std::endl;

#else
    // save/visualize initial solution
    model_solver.visualize(filename, 0, 0., subsampling);

    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    // auto stokes_vec = model_solver.stokes_vec();
    // auto pfield_vec = model_solver.pfield_vec(0);
    // auto ofield_vec = model_solver.ofield_vec(0);


    const size_t num_cells = model_solver.num_cells_;
    std::vector<std::vector<typename CellModelSolver::VectorType>> ret(1 + 2 * num_cells);
    for (size_t kk = 0; kk < num_cells; ++kk) {
      ret[kk].push_back(model_solver.pfield_vec(kk));
      ret[num_cells + kk].push_back(model_solver.ofield_vec(kk));
    }
    ret[2 * num_cells].push_back(model_solver.stokes_vec());

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t << std::endl;
      for (size_t kk = 0; kk < num_cells; ++kk) {
        // model_solver.set_pfield_vec(kk, pfield_vec);
        // model_solver.set_ofield_vec(kk, ofield_vec);
        // model_solver.set_stokes_vec(stokes_vec);
        // model_solver.prepare_pfield_op(dt, kk);
        // pfield_vec = model_solver.apply_inverse_pfield_op_with_param(pfield_vec, kk, Be, Ca, Pa);
        // std::cout << "Pfield " << kk << " done" << std::endl;
        // model_solver.set_pfield_vec(kk, pfield_vec);
        // model_solver.set_ofield_vec(kk, ofield_vec);
        // model_solver.set_stokes_vec(stokes_vec);
        // model_solver.prepare_ofield_op(dt, kk);
        // ofield_vec = model_solver.apply_inverse_ofield_op_with_param(ofield_vec, kk, Pa);
        // model_solver.set_pfield_vec(kk, pfield_vec);
        // model_solver.set_ofield_vec(kk, ofield_vec);
        // std::cout << "Ofield " << kk << " done" << std::endl;
        model_solver.set_pfield_vec(kk, ret[kk].back());
        model_solver.prepare_pfield_op(dt, kk);
        ret[kk].push_back(model_solver.apply_inverse_pfield_op_with_param(ret[kk].back(), kk, Be, Ca, Pa));
        std::cout << "Pfield " << kk << " done" << std::endl;
        model_solver.set_pfield_vec(kk, ret[kk].back());
        model_solver.prepare_ofield_op(dt, kk);
        ret[num_cells + kk].push_back(
            model_solver.apply_inverse_ofield_op_with_param(ret[kk + num_cells].back(), kk, Pa));
        model_solver.set_ofield_vec(kk, ret[num_cells + kk].back());
        std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      model_solver.prepare_stokes_op();
      ret[2 * num_cells].push_back(model_solver.apply_inverse_stokes_op());
      model_solver.set_stokes_vec(ret[2 * num_cells].back());
      std::cout << "Stokes done" << std::endl;

      t += actual_dt;

      // check if data should be written in this timestep (and write)
      //  if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t, next_save_time)) {
      //    model_solver.visualize(filename, save_step_counter, t, subsampling);
      //    next_save_time += write_step;
      //    ++save_step_counter;
      //  }
    } // while (t < t_end)
    for (size_t ii = 0; ii < ret[0].size(); ++ii) {
      model_solver.set_pfield_vec(0, ret[0][ii]);
      model_solver.set_ofield_vec(0, ret[num_cells][ii]);
      model_solver.set_stokes_vec(ret[2 * num_cells][ii]);
      model_solver.visualize(filename, ii, ii, subsampling);
    }
#endif
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
