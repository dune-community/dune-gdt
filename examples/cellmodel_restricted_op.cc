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

#include <random>

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
                                 false,
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
    CellModelSolver model_solver2(testcase,
                                  t_end,
                                  num_elements_x,
                                  num_elements_y,
                                  false,
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

    const size_t max_output_dof = model_solver.pfield_vector().size();
    const size_t num_output_dofs = 20;
    std::vector<size_t> pfield_output_dofs(num_output_dofs);
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<typename std::mt19937::result_type> distrib(0, max_output_dof);
    for (size_t ii = 0; ii < num_output_dofs; ++ii)
      pfield_output_dofs[ii] = distrib(rng);
    const size_t num_cells = model_solver.num_cells_;
    auto source = model_solver.pfield_vector();
    for (size_t kk = 0; kk < num_cells; ++kk) {
      model_solver.prepare_restricted_pfield_operator(pfield_output_dofs, dt, kk);
      const auto& pfield_input_dofs = model_solver.pfield_restricted_op_input_dofs(kk);
      const size_t num_input_dofs = pfield_input_dofs.size();
      typename CellModelSolver::VectorType restricted_source(num_input_dofs, 0.);
      for (size_t ii = 0; ii < num_input_dofs; ++ii)
        restricted_source[ii] = source[pfield_input_dofs[ii]];
      auto restricted_result = model_solver.apply_restricted_pfield_operator(restricted_source, kk);
      model_solver2.prepare_pfield_operator(dt, kk);
      auto result = model_solver2.apply_pfield_operator(source, kk);
      for (size_t ii = 0; ii < num_output_dofs; ++ii)
        if (XT::Common::FloatCmp::ne(restricted_result[ii], result[pfield_output_dofs[ii]], 1e-13, 1e-13))
          std::cout << "Failed: " << ii << ", " << pfield_output_dofs[ii] << ", " << result[pfield_output_dofs[ii]]
                    << ", " << restricted_result[ii] << std::endl;
    }
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
