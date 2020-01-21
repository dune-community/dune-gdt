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

#include <chrono>
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

    // GMRES options
    const double gmres_reduction = DXTC_CONFIG_GET("gmres_reduction", 1e-10);
    const int direct_gmres_restart = DXTC_CONFIG_GET("direct_gmres_restart", 100);
    const int custom_gmres_restart = DXTC_CONFIG_GET("custom_gmres_restart", 40);
    const int gmres_verbose = DXTC_CONFIG_GET("gmres_verbose", 0);

    CellModelSolver model_solver(testcase,
                                 t_end,
                                 num_elements_x,
                                 num_elements_y,
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
                                 "custom",
                                 "schur",
                                 gmres_reduction,
                                 custom_gmres_restart,
                                 gmres_verbose,
                                 linearize);
    CellModelSolver model_solver2(testcase,
                                  t_end,
                                  num_elements_x,
                                  num_elements_y,
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
                                  "direct",
                                  "schur",
                                  gmres_reduction,
                                  direct_gmres_restart,
                                  gmres_verbose,
                                  linearize);

    // CellModelSolver model_solver3(testcase,
    //                               t_end,
    //                               num_elements_x,
    //                               num_elements_y,
    //                               false,
    //                               Be,
    //                               Ca,
    //                               Pa,
    //                               Re,
    //                               Fa,
    //                               xi,
    //                               kappa,
    //                               c_1,
    //                               beta,
    //                               gamma,
    //                               epsilon,
    //                               In,
    //                               "custom",
    //                               "direct",
    //                               linearize);


    std::chrono::duration<double> ref_time(0.);
    std::chrono::duration<double> pfield_direct_time(0.);
    std::chrono::duration<double> ofield_direct_time(0.);
    auto begin = std::chrono::steady_clock::now();
    const auto result1 = model_solver.solve(dt, false, write_step, filename, subsampling);
    ref_time = std::chrono::steady_clock::now() - begin;
    begin = std::chrono::steady_clock::now();
    const auto result2 = model_solver2.solve(dt, false, write_step, filename, subsampling);
    pfield_direct_time = std::chrono::steady_clock::now() - begin;
    // begin = std::chrono::steady_clock::now();
    // const auto result3 = model_solver3.solve(dt, false, write_step, filename, subsampling);
    // ofield_direct_time = std::chrono::steady_clock::now() - begin;
    std::cout << "Timings: ref =  " << ref_time.count() << "  vs. pfield direct = " << pfield_direct_time.count()
              << " vs. ofield direct " << ofield_direct_time.count() << std::endl;
    for (size_t ii = 0; ii < result1.size(); ++ii) {
      for (size_t jj = 0; jj < result1[ii].size(); ++jj) {
        for (size_t kk = 0; kk < result1[ii][jj].size(); ++kk) {
          if (XT::Common::FloatCmp::ne(result1[ii][jj][kk], result2[ii][jj][kk], 1e-8, 1e-8))
            std::cout << "Entries (" << ii << ", " << jj << ", " << kk << ") with different pfield solvers differ by "
                      << std::abs(result1[ii][jj][kk] - result2[ii][jj][kk]) << ", " << result1[ii][jj][kk] << " vs. "
                      << result1[ii][jj][kk] << std::endl;
          // if (XT::Common::FloatCmp::ne(result1[ii][jj][kk], result3[ii][jj][kk], 1e-8, 1e-8))
          // std::cout << "Entries (" << ii << ", " << jj << ", " << kk << ") with different ofield solvers differ by "
          // << std::abs(result1[ii][jj][kk] - result2[ii][jj][kk]) << ", " << result1[ii][jj][kk] << " vs. "
          // << result1[ii][jj][kk] << std::endl;
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
