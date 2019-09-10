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
    std::string filename = config.get("output.filename", "drosophila") + (linearize ? "_linearized" : "");
    bool subsampling = config.get<bool>("output.subsampling", true);
    // a negative value of write step is interpreted as "write all steps"
    double write_step = config.template get<double>("output.write_step", -1.);

    // create grid for [0, 160] x [0, 40] with periodic boundaries in x-direction
    FieldVector<double, d> lower_left, upper_right;
    std::string periodic_dirs;
    size_t num_cells;
    if (testcase == "single_cell") {
      lower_left = {{0., 0.}};
      upper_right = {{160., 40.}};
      periodic_dirs = "01";
      num_cells = 1;
    } else if (testcase == "two_cells") {
      lower_left = {{0., 0.}};
      upper_right = {{50., 50.}};
      // periodic_dirs = "11";
      periodic_dirs = "00";
      num_cells = 2;
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    }
    auto grid = XT::Grid::make_cube_grid<G>(lower_left, upper_right, /*num_elements=*/{num_elements_x, num_elements_y});
    grid.grid().globalRefine(1);
    auto nonperiodic_grid_view = grid.leaf_view();
    std::bitset<d> periodic_directions(periodic_dirs);
    PGV grid_view(nonperiodic_grid_view, periodic_directions);
    std::shared_ptr<const XT::Functions::FunctionInterface<d, d>> u_initial_func;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> phi_initial_funcs;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> mu_initial_funcs;
    std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d, d>>> P_initial_funcs;

    if (testcase == "single_cell") {
      // create and project initial values
      // we only need initial values for P and phi
      // mu_initial is only needed if linearization is used
      // Initially, cell is circular with Radius R=5 and placed in the center of the domain
      // \Omega = [0, 160] \times [0, 40].
      // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
      // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
      FieldVector<double, d> center{upper_right[0] / 2., upper_right[1] / 2.};
      auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
      phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
          50,
          /*evaluate=*/
          [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2.) * epsilon)); },
          /*name=*/"phi_initial"));
      mu_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d>>(
          50,
          /*evaluate=*/
          [& phi_in = phi_initial_funcs[0], epsilon](const auto& x, const auto& param) {
            // TODO: add approximation of laplacian term
            const auto phi = phi_in->evaluate(x, param);
            return 1. / epsilon * (std::pow(phi, 3) - phi);
          },
          /*name=*/"mu_initial"));

      // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
      // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
      std::srand(1); // set seed for std::rand to 1
      P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
          50,
          /*evaluate=*/
          [& phi_in = phi_initial_funcs[0]](const auto& x, const auto& param) {
            // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
            // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
            // auto ret = FieldVector<double, d>({1. + rand1, 0. +
            // rand2});
            auto ret = FieldVector<double, d>({1., 0.});
            ret *= (phi_in->evaluate(x, param) + 1.) / 2.;
            return ret;
          },
          /*name=*/"P_initial"));
      u_initial_func = std::make_shared<const XT::Functions::ConstantFunction<d, d>>(0.);

      // interpolate initial and boundary values
    } else if (testcase == "two_cells") {
      // create and project initial values
      // we only need initial values for P_i and phi_i
      // mu_initial is only needed if linearization is used
      // Initially, cell is circular with Radius R=5 and placed in the center of the domain
      // \Omega = [0, 160] \times [0, 40].
      // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
      // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
      FieldVector<double, d> center1{15, 15};
      FieldVector<double, d> center2{35, 35};
      auto r1 = [center1](const auto& xr) { return 4.0 - (center1 - xr).two_norm(); };
      auto r2 = [center2](const auto& xr) { return 4.0 - (center2 - xr).two_norm(); };
      const XT::Functions::GenericFunction<d> phi1_initial(50,
                                                           /*evaluate=*/
                                                           [r = r1, epsilon](const auto& x, const auto& /*param*/) {
                                                             return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
                                                           },
                                                           /*name=*/"phi1_initial");
      const XT::Functions::GenericFunction<d> mu1_initial(50,
                                                          /*evaluate=*/
                                                          [phi1_initial, epsilon](const auto& x, const auto& param) {
                                                            // TODO: add approximation of laplacian term
                                                            const auto phi = phi1_initial.evaluate(x, param);
                                                            return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                          },
                                                          /*name=*/"mu1_initial");
      const XT::Functions::GenericFunction<d> phi2_initial(50,
                                                           /*evaluate=*/
                                                           [r = r2, epsilon](const auto& x, const auto& /*param*/) {
                                                             return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
                                                           },
                                                           /*name=*/"phi1_initial");
      const XT::Functions::GenericFunction<d> mu2_initial(50,
                                                          /*evaluate=*/
                                                          [phi2_initial, epsilon](const auto& x, const auto& param) {
                                                            // TODO: add approximation of laplacian term
                                                            const auto phi = phi2_initial.evaluate(x, param);
                                                            return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                          },
                                                          /*name=*/"mu1_initial");

      // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
      // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
      std::srand(1); // set seed for std::rand to 1
      const XT::Functions::GenericFunction<d, d> P1_initial(50,
                                                            /*evaluate=*/
                                                            [phi1_initial](const auto& x, const auto& param) {
                                                              // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto ret = FieldVector<double, d>({1. + rand1, 0. +
                                                              // rand2});
                                                              auto ret = FieldVector<double, d>({1., 0.});
                                                              ret *= (phi1_initial.evaluate(x, param) + 1.) / 2.;
                                                              return ret;
                                                            },
                                                            /*name=*/"P_initial");
      const XT::Functions::GenericFunction<d, d> P2_initial(50,
                                                            /*evaluate=*/
                                                            [phi2_initial](const auto& x, const auto& param) {
                                                              // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                              // auto ret = FieldVector<double, d>({1. + rand1, 0. +
                                                              // rand2});
                                                              auto ret = FieldVector<double, d>({1., 0.});
                                                              ret *= (phi2_initial.evaluate(x, param) + 1.) / 2.;
                                                              return ret;
                                                            },
                                                            /*name=*/"P_initial");

      const XT::Functions::ConstantFunction<d, d> u_initial(0.);
    } else {
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    }

    // interpolate initial and boundary values
    auto phi_space = make_continuous_lagrange_space<1>(grid_view, 2);
    auto u_space = make_continuous_lagrange_space<d>(grid_view, 2);
    VectorType u_initial_vec(u_space.mapper().size());
    std::vector<VectorType> phi_initial_vecs(num_cells, VectorType(phi_space.mapper().size()));
    std::vector<VectorType> mu_initial_vecs(num_cells, VectorType(phi_space.mapper().size()));
    std::vector<VectorType> P_initial_vecs(num_cells, VectorType(u_space.mapper().size()));
    auto u_initial = make_discrete_function(u_space, u_initial_vec);
    const XT::Functions::ConstantFunction<d> minus_one(-1.);
    default_interpolation(*u_initial_func, u_initial);
    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
    for (size_t kk = 0; kk < num_cells; kk++) {
      auto phi_initial = make_discrete_function(phi_space, phi_initial_vecs[kk]);
      default_interpolation(*phi_initial_funcs[kk], phi_initial);
      boundary_interpolation(minus_one, phi_initial, all_dirichlet_boundary_info, XT::Grid::DirichletBoundary{});
      auto mu_initial = make_discrete_function(phi_space, mu_initial_vecs[kk]);
      default_interpolation(*mu_initial_funcs[kk], mu_initial);
      auto P_initial = make_discrete_function(u_space, P_initial_vecs[kk]);
      default_interpolation(*P_initial_funcs[kk], P_initial);
    }

    static const size_t pol_order = 2;
    CellModelSolver model_solver(u_initial_vec,
                                 P_initial_vecs,
                                 phi_initial_vecs,
                                 mu_initial_vecs,
                                 num_cells,
                                 lower_left,
                                 upper_right,
                                 periodic_dirs,
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
                                 linearize,
                                 pol_order);

    // save/visualize initial solution
    model_solver.visualize(filename, 0, 0., subsampling);

    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    VectorType& stokes_vec = model_solver.stokes_vector();
    VectorType& pfield_vec = model_solver.pfield_vector();
    VectorType& ofield_vec = model_solver.ofield_vector();

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
        pfield_vec = model_solver.solve_pfield_system(pfield_vec, kk);
        model_solver.set_phasefield_variables(kk, pfield_vec);
        std::cout << "Pfield " << kk << " done" << std::endl;
        model_solver.prepare_ofield_operator(dt, kk);
        ofield_vec = model_solver.solve_ofield_system(ofield_vec, kk);
        model_solver.set_orientationfield_variables(kk, ofield_vec);
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
