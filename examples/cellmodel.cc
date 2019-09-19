// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include "cellmodel.hh"

int main(int argc, char* argv[])
{
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
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;
    std::cout << "Ca: " << Ca << ", Be: " << Be << ", Pa: " << Pa << ", Fa: " << Fa << ", Re: " << Re << std::endl;

    // output
    std::string filename = config.get("output.filename", "cellmodel") + (linearize ? "_linearized" : "");
    bool subsampling = config.get<bool>("output.subsampling", true);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);

    CellModelSolver model_solver(testcase,
                                 t_end,
                                 num_elements_x,
                                 num_elements_y,
                                 Re,
                                 Fa,
                                 xi,
                                 kappa,
                                 c_1,
                                 Pa,
                                 beta,
                                 gamma,
                                 Be,
                                 Ca,
                                 epsilon,
                                 In,
                                 linearize);
    auto result = model_solver.solve(dt, write_step, filename, subsampling);

    for (const auto& vec1 : result[0])
      std::cout << &vec1 << std::endl;

#if 0
    // save/visualize initial solution
    model_solver.visualize(filename, 0, 0., subsampling);

    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    VectorType stokes_vec = model_solver.stokes_vector();
    VectorType pfield_vec = model_solver.pfield_vector();
    VectorType ofield_vec = model_solver.ofield_vector();

    const size_t num_cells = model_solver.num_cells_;
    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      std::cout << "Current time: " << t << std::endl;
      for (size_t kk = 0; kk < num_cells; ++kk) {
        model_solver.prepare_pfield_operator(dt, kk);
        pfield_vec = model_solver.solve_pfield(pfield_vec, kk);
        model_solver.set_pfield_variables(kk, pfield_vec);
        std::cout << "Pfield " << kk << " done" << std::endl;
        model_solver.prepare_ofield_operator(dt, kk);
        ofield_vec = model_solver.solve_ofield(ofield_vec, kk);
        model_solver.set_ofield_variables(kk, ofield_vec);
        std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      model_solver.prepare_stokes_operator();
      stokes_vec = model_solver.solve_stokes();
      model_solver.set_stokes_variables(stokes_vec);
      std::cout << "Stokes done" << std::endl;

      t += actual_dt;

      // check if data should be written in this timestep (and write)
      if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t, next_save_time)) {
        model_solver.visualize(filename, save_step_counter, t, subsampling);
        next_save_time += write_step;
        ++save_step_counter;
      }
    } // while (t < t_end)
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
