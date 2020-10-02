// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014, 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2016)

#include "config.h"

#include "boltzmann.hh"

int main(int argc, char* argv[])
{
  static const size_t dimDomain = 2;
  static const bool linear = true;
  try {
    // parse options
    if (argc == 1)
      std::cout << "The following options are available: " << argv[0]
                << " [-output_dir DIR -num_save_steps INT -gridsize INT "
                << "  -sigma_s_1 FLOAT -sigma_s_2 FLOAT -sigma_a_1 FLOAT -sigma_a_2 FLOAT"
                << " --no_visualization --silent --random_parameters]" << std::endl;

    size_t num_save_steps = 10;
    size_t grid_size = 20;
    bool visualize = true;
    bool silent = false;
    bool random_parameters = false;
    bool parameters_given = false;
    std::string output_dir;
    double sigma_s_lower = 0, sigma_s_upper = 8, sigma_a_lower = 0, sigma_a_upper = 8;
    double sigma_s_1 = 1, sigma_s_2 = 0, sigma_a_1 = 0, sigma_a_2 = 10;
    for (int i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "-output_dir") {
        if (i + 1 < argc) {
          output_dir = argv[++i];
        } else {
          std::cerr << "-output_dir option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-num_save_steps") {
        if (i + 1 < argc) {
          num_save_steps = XT::Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-num_save_steps option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "--no_visualization") {
        visualize = false;
      } else if (std::string(argv[i]) == "--silent") {
        silent = true;
      } else if (std::string(argv[i]) == "--random_parameters") {
        if (parameters_given) {
          std::cerr << "You specified a value for at least one parameter so you can't use --random_parameters!"
                    << std::endl;
          return 1;
        }
        random_parameters = true;
        RandomNumberGeneratorType rng{std::random_device()()};
        std::uniform_real_distribution<double> sigma_s_dist(sigma_s_lower, sigma_s_upper);
        std::uniform_real_distribution<double> sigma_a_dist(sigma_a_lower, sigma_a_upper);
        sigma_s_1 = sigma_s_dist(rng);
        sigma_s_2 = sigma_s_dist(rng);
        sigma_a_1 = sigma_a_dist(rng);
        sigma_a_2 = sigma_a_dist(rng);
      } else if (std::string(argv[i]) == "-gridsize") {
        if (i + 1 < argc) {
          grid_size = XT::Common::from_string<size_t>(argv[++i]);
        } else {
          std::cerr << "-gridsize option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_1 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_s_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_s_2 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_s_2 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_1") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_1 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_1 option requires one argument." << std::endl;
          return 1;
        }
      } else if (std::string(argv[i]) == "-sigma_a_2") {
        if (random_parameters) {
          std::cerr << "You specified a value for at least one parameter on the command line so you can't use "
                    << "--random_parameters!" << std::endl;
          return 1;
        }
        if (i + 1 < argc) {
          sigma_a_2 = XT::Common::from_string<double>(argv[++i]);
          parameters_given = true;
        } else {
          std::cerr << "-sigma_a_2 option requires one argument." << std::endl;
          return 1;
        }
      } else {
        std::cerr << "Unknown option " << std::string(argv[i]) << std::endl;
        return 1;
      }
    }

    std::ofstream parameterfile;
    parameterfile.open(output_dir + "_parameters.txt");
    parameterfile << "Gridsize: " << XT::Common::to_string(grid_size) + " x " + XT::Common::to_string(grid_size)
                  << std::endl;

    // run solver
    parameterfile << "Domain was composed of two materials, parameters were: " << std::endl
                  << "First material: sigma_s = " + XT::Common::to_string(sigma_s_1)
                         + ", sigma_a = " + XT::Common::to_string(sigma_a_1)
                  << std::endl
                  << "Second material: sigma_s = " + XT::Common::to_string(sigma_s_2)
                         + ", sigma_a = " + XT::Common::to_string(sigma_a_2)
                  << std::endl;

    auto solver = std::make_shared<BoltzmannSolver<dimDomain>>(
        output_dir, num_save_steps, grid_size, visualize, silent, sigma_s_1, sigma_s_2, sigma_a_1, sigma_a_2);

    DXTC_TIMINGS.start("solve_all");
    if (!linear) {
      std::vector<size_t> output_dofs{2728, 3868, 4468, 929};
      solver->prepare_restricted_operator(output_dofs);
      using VectorType = typename XT::LA::Container<double, XT::LA::Backends::common_dense>::VectorType;
      auto initial_vals = solver->get_initial_values();
      RandomNumberGeneratorType rng{std::random_device()()};
      std::uniform_real_distribution<double> distribution(1, 1000);
      for (auto&& val : initial_vals)
        val *= 1e4 * distribution(rng);
      auto source_dofs = solver->restricted_op_input_dofs();
      std::cout << source_dofs << std::endl;
#if 0
    VectorType initial_vals_restr{
        3.00887845e-05, 7.40090567e-05, 7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 1.00000000e-08,
        3.79301179e-05, 1.42264780e-05, 1.51590332e-05, 1.00000000e-08, 6.04617301e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 3.00887845e-05, 1.00000000e-08, 3.39443780e-04, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05,
        7.40090567e-05, 1.00000000e-08, 3.39443780e-04, 1.51590332e-05, 1.42264780e-05, 3.79301179e-05, 1.00000000e-08,
        1.00000000e-08, 6.04617301e-05, 1.47623908e-05, 9.55564780e-06, 9.55564780e-06, 1.47623908e-05, 1.00000000e-08,
        3.30163508e-05, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.20284294e-04, 1.00000000e-08, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08,
        1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02,
        1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08, 1.20284294e-04,
        1.00000000e-08, 1.00000000e-08, 3.62665787e-02, 3.68187108e-02, 3.62665787e-02, 3.68187108e-02, 3.68187108e-02,
        3.62665787e-02, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.50692562e-04, 2.63568632e-05, 6.74968960e-05, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 2.63568632e-05, 6.74968960e-05, 1.00000000e-08,
        1.50692562e-04, 2.21867566e-04, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.20284294e-04, 1.00000000e-08, 1.00000000e-08,
        1.00000000e-08, 3.00887845e-05, 7.40090567e-05, 3.00887845e-05, 7.40090567e-05, 3.39443780e-04, 1.00000000e-08};
#endif
      VectorType initial_vals_restr(source_dofs.size());
      for (size_t kk = 0; kk < source_dofs.size(); ++kk)
        initial_vals_restr[kk] = initial_vals[source_dofs[kk]];
      auto output1 = solver->apply_restricted_kinetic_operator(initial_vals_restr);
      // for (size_t ii = 0; ii < 1000; ++ii) {
      const size_t ii = 0;
      const auto actual_initial_vals = initial_vals_restr * ii / 1000.;
      output1 = solver->apply_restricted_kinetic_operator(actual_initial_vals);
      // }
      auto output2 = solver->apply_kinetic_operator(initial_vals, 0, solver->time_step_length());
      for (size_t kk = 0; kk < output_dofs.size(); ++kk)
        EXPECT_NEAR(output1[kk], output2[output_dofs[kk]], 1e-10);
    }
    const auto result = solver->solve();
    std::cout << " Result = " << std::accumulate(result.back().begin(), result.back().end(), 0.) << std::endl;
    DXTC_TIMINGS.stop("solve_all");
    parameterfile << "Elapsed time: " << DXTC_TIMINGS.walltime("solve_all") / 1000.0 << " s" << std::endl;
    parameterfile.close();

    return 0;
  } catch (Dune::Exception& e) {
    std::cerr << "Dune reported: " << e.what() << std::endl;
    std::abort();
  }
} // ... main(...)