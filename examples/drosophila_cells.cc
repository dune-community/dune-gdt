// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/solver/istl/saddlepoint.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>

#include <dune/gdt/interpolations/default.hh>


using namespace Dune;
using namespace Dune::GDT;


// some global defines
using G = YASP_2D_EQUIDISTANT_OFFSET;
static const constexpr size_t d = G::dimension;
using GV = typename G::LeafGridView;
using PGV = XT::Grid::PeriodicGridView<GV>;
using E = XT::Grid::extract_entity_t<GV>;
using I = XT::Grid::extract_intersection_t<GV>;
using PI = XT::Grid::extract_intersection_t<PGV>;
using MatrixType = XT::LA::IstlRowMajorSparseMatrix<double>;
using VectorType = XT::LA::IstlDenseVector<double>;

template <class DiscreteFunc>
static std::string get_filename(const std::string& prefix, const DiscreteFunc& func, const size_t step)
{
  return prefix + "_" + func.name() + XT::Common::to_string(step) + ".txt";
}

template <class DiscreteFunc>
static std::string
get_rankfilename(const std::string& prefix, const DiscreteFunc& func, const size_t step, const int rank)
{
  return prefix + "_rank_" + XT::Common::to_string(rank) + "_" + func.name() + XT::Common::to_string(step) + ".txt";
}

template <class DiscreteFunc>
static void write_to_textfile(const DiscreteFunc& func, const std::string& prefix, const size_t step, const double t)
{
  const auto& grid_view = func.space().grid_view();
  const auto local_func = func.local_function();
  // write one file per MPI rank
  std::ofstream rankfile(get_rankfilename(prefix, func, step, grid_view.comm().rank()));
  for (const auto& entity : elements(grid_view, Dune::Partitions::interiorBorder)) {
    local_func->bind(entity);
    const auto entity_center = entity.geometry().center();
    auto position = entity_center;
    assert(position.size() == d);
    const auto val = local_func->evaluate(entity.geometry().local(position), {"t", t});
    for (size_t ii = 0; ii < d; ++ii)
      rankfile << XT::Common::to_string(position[ii], 15) << " ";
    rankfile << val << std::endl;
  }
  rankfile.close();
  // Wait till files on all ranks are written
  grid_view.comm().barrier();
  // Merge files
  if (grid_view.comm().rank() == 0) {
    const std::string merged_file_name = get_filename(prefix, func, step);
    std::remove(merged_file_name.c_str());
    std::ofstream merged_file(merged_file_name, std::ios_base::binary | std::ios_base::app);
    for (int ii = 0; ii < grid_view.comm().size(); ++ii) {
      const std::string rankfile_to_merge_name = get_rankfilename(prefix, func, step, ii);
      std::ifstream rankfile_to_merge(rankfile_to_merge_name, std::ios_base::binary);
      merged_file << rankfile_to_merge.rdbuf();
      rankfile_to_merge.close();
      std::remove(rankfile_to_merge_name.c_str());
    } // ii
    merged_file.close();
  } // if (rank == 0)
} // void write_to_textfile()

template <class DiscreteFunc, class VecDiscreteFunc>
static void write_files(const bool visualize,
                        const VecDiscreteFunc& u,
                        const DiscreteFunc& p,
                        const VecDiscreteFunc& P,
                        const VecDiscreteFunc& Pnat,
                        const DiscreteFunc& phi,
                        const DiscreteFunc& phinat,
                        const DiscreteFunc& mu,
                        const std::string& prefix,
                        const size_t step,
                        const double t)
{
  std::string postfix = "_" + Dune::XT::Common::to_string(step);
  if (visualize) {
    u.visualize(prefix + "_" + u.name() + postfix);
    p.visualize(prefix + "_" + p.name() + postfix);
    P.visualize(prefix + "_" + P.name() + postfix);
    Pnat.visualize(prefix + "_" + Pnat.name() + postfix);
    phi.visualize(prefix + "_" + phi.name() + postfix);
    phinat.visualize(prefix + "_" + phinat.name() + postfix);
    mu.visualize(prefix + "_" + mu.name() + postfix);
  }
  write_to_textfile(u, prefix, step, t);
  write_to_textfile(p, prefix, step, t);
  write_to_textfile(P, prefix, step, t);
  write_to_textfile(Pnat, prefix, step, t);
  write_to_textfile(phi, prefix, step, t);
  write_to_textfile(phinat, prefix, step, t);
  write_to_textfile(mu, prefix, step, t);
}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_navier_stokes(VecDiscreteFunc& u,
                         DiscreteFunc& p,
                         const VecDiscreteFunc& P,
                         const VecDiscreteFunc& Pnat,
                         const DiscreteFunc& phi,
                         const DiscreteFunc& phinat,
                         const DiscreteFunc& mu,
                         const double Re,
                         const double Fa,
                         const double xi,
                         const XT::Grid::BoundaryInfo<PI>& boundary_info)
{
  const auto& u_space = u.space();
  const auto& p_space = p.space();

  const auto& grid_view = u_space.grid_view();
  // Setup spaces and matrices and vectors
  // Equations are
  // \int \nabla u \nabla v - \int p div v = \int ff v
  // \int (div u) q = \int gg q
  // System is [A B; B^T C] [u; p] = [f; g]
  // Dimensions are: A: n x n, B: n x m, C: m x m, u: n, f: n, p: m, g: m
  const size_t m = u_space.mapper().size();
  const size_t n = p_space.mapper().size();
  auto pattern_A = make_element_sparsity_pattern(u_space, u_space, grid_view);
  auto pattern_B = make_element_sparsity_pattern(u_space, p_space, grid_view);
  auto pattern_C = make_element_sparsity_pattern(p_space, p_space, grid_view);
  MatrixType A(m, m, pattern_A);
  MatrixType B(m, n, pattern_B);
  MatrixOperator<MatrixType, PGV, d> A_operator(grid_view, u_space, u_space, A);
  MatrixOperator<MatrixType, PGV, 1, 1, d> B_operator(grid_view, p_space, u_space, B);
  // calculate A_{ij} as \int \nabla v_i \nabla v_j
  A_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalEllipticIntegrand<E, d>(1.)));
  // calculate B_{ij} as \int \nabla p_i div(v_j)
  B_operator.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
      LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
  // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
  VectorType f_vector(m), g_vector(n, 0.), p_basis_integrated_vector(n);
  auto f_functional = make_vector_functional(u_space, f_vector);
  // TODO: implement rhs
  XT::Functions::ConstantGridFunction<E, d> f(0.);
  f_functional.append(LocalElementIntegralFunctional<E, d>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), f)));
  A_operator.append(f_functional);
  auto p_basis_integrated_functional = make_vector_functional(p_space, p_basis_integrated_vector);
  XT::Functions::ConstantGridFunction<E> one_function(1);
  p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), one_function)));
  B_operator.append(p_basis_integrated_functional);
  // Dirichlet constrainst for u
  auto dirichlet_constraints = make_dirichlet_constraints(u_space, boundary_info);
  A_operator.append(dirichlet_constraints);
  // assemble everything
  A_operator.assemble(true);
  B_operator.assemble(true);
  // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
  // Dirichlet DoF, so we also need to clear the row of B.
  dirichlet_constraints.apply(A, f_vector);
  for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs())
    B.clear_row(DoF);

  // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
  // [A B; B^T 0] to the unit vector.
  MatrixType C(n, n, pattern_C);
  size_t dof_index = 0;
  B.clear_col(dof_index);
  g_vector.set_entry(dof_index, 0.);
  C.set_entry(dof_index, dof_index, 1.);

  // now solve the system
  XT::LA::SaddlePointSolver<double> solver(A, B, B, C);
  // solve both by direct solver and by schurcomplement (where the schur complement is inverted by CG and the inner
  // solves with A are using a direct method)
  std::string type = "cg_direct_schurcomplement";
  solver.apply(f_vector, g_vector, u.dofs().vector(), p.dofs().vector(), type);

  // ensure int_\Omega p = 0
  auto p_integral = p_basis_integrated_vector * p.dofs().vector();
  auto p_correction = make_discrete_function<VectorType>(p_space, "p_corr");
  auto vol_domain = 4.;
  XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain);
  interpolate(const_p_integral_func, p_correction);
  p -= p_correction;
}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_ofield(const VecDiscreteFunc& u,
                  const DiscreteFunc& p,
                  VecDiscreteFunc& P,
                  VecDiscreteFunc& Pnat,
                  const DiscreteFunc& phi,
                  const DiscreteFunc& phinat,
                  const DiscreteFunc& mu,
                  const double xi,
                  const double kappa,
                  const double c_1,
                  const double Pa,
                  const double beta)
{}

template <class DiscreteFunc, class VecDiscreteFunc>
void solve_pfield(const VecDiscreteFunc& u,
                  const DiscreteFunc& p,
                  const VecDiscreteFunc& P,
                  const VecDiscreteFunc& Pnat,
                  DiscreteFunc& phi,
                  DiscreteFunc& phinat,
                  DiscreteFunc& mu,
                  const double gamma,
                  const double c_1,
                  const double Pa,
                  const double beta,
                  const double epsilon)
{}

int main(int argc, char* argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    if (argc > 1)
      DXTC_CONFIG.read_options(argc, argv);
    XT::Common::TimedLogger().create(DXTC_CONFIG_GET("logger.info", 1), DXTC_CONFIG_GET("logger.debug", -1));
    auto logger = XT::Common::TimedLogger().get("main");

    // read configuration
    XT::Common::Configuration config("activepolargels.ini");

    // grid config
    unsigned int num_elements_x = config.template get<unsigned int>("grid.NX", static_cast<unsigned int>(16));
    unsigned int num_elements_y = config.template get<unsigned int>("grid.NY", static_cast<unsigned int>(4));
    std::string periodic_dirs = config.get("grid.periodic_dirs", "01");

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
    double Re = rho * U * L / eta;
    double Ca = 2. * std::sqrt(2) / 3. * eta * U / sigma;
    double Be = 4. * std::sqrt(2) / 3. * eta * U * L * L / b_N;
    double Pa = eta * U * L / k;
    double Fa = eta * U / (zeta * L);
    const double kappa = eta_rot / eta;

    // output
    std::string filename = config.get("output.filename", "drosophila");
    double write_step = config.template get<double>("output.write_step", 0.01);

    // create grid for [0, 160] x [0, 40] with periodic boundaries in x-direction
    FieldVector<double, d> upper_right{160., 40.};
    auto grid = XT::Grid::make_cube_grid<G>(
        /*lower_left=*/{0., 0.}, upper_right, /*num_elements=*/{num_elements_x, num_elements_y});
    auto nonperiodic_grid_view = grid.leaf_view();
    std::bitset<d> periodic_directions(periodic_dirs);
    PGV grid_view(nonperiodic_grid_view, periodic_directions);

    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    XT::Grid::AllDirichletBoundaryInfo<PI> dirichlet_boundary_info;
    XT::Grid::AllNeumannBoundaryInfo<PI> neumann_boundary_info;

    const XT::Functions::ConstantFunction<d> lambda(1);

    // create spaces for the variables
    // Taylor-Hood P2-P1 space for the Stokes variables (u, p) \in S_{stokes} = (H^1_0)^d x L^2\R, test functions (v, q)
    // \in S_{stokes}
    auto u_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    auto p_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/1);

    // P2 spaces for the orientation field variables P, P^{\natural} \in S_{ofield} = (H^1)^d x (L^2)^d, test functions
    // (Q, Q^\natural) \in S_{ofield}
    auto P_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    auto Pnat_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);

    // P2 spaces for the phase field variables \phi, \phi^\natural, \mu \in S_{pfield} = H^1_0 x H^1 x H^1, test
    // functions (\chi, \chi^\natural, \nu) \in S_{pfield}
    auto phi_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);
    auto phinat_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);
    auto mu_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);

    // create discrete functions for the variables
    auto u = make_discrete_function<VectorType>(u_space, "u");
    auto p = make_discrete_function<VectorType>(p_space, "p");
    auto P = make_discrete_function<VectorType>(P_space, "P");
    auto Pnat = make_discrete_function<VectorType>(Pnat_space, "Pnat");
    auto phi = make_discrete_function<VectorType>(phi_space, "phi");
    auto phinat = make_discrete_function<VectorType>(phinat_space, "phinat");
    auto mu = make_discrete_function<VectorType>(mu_space, "mu");

    // create and project initial values
    // we only need initial values for P and phi
    // Initially, cell is circular with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 160] \times [0, 40].
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
    FieldVector<double, d> center{upper_right[0] / 2, upper_right[1] / 2};
    auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
    const XT::Functions::GenericFunction<d> phi_initial(
        10,
        /*evaluate=*/
        [r, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2) * epsilon)); },
        /*name=*/"phi_initial");

    // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
    // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
    std::srand(1); // set seed for std::rand to 1
    const XT::Functions::GenericFunction<d, d> P_initial(3,
                                                         /*evaluate=*/
                                                         [phi_initial](const auto& x, const auto& param) {
                                                           auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
                                                           auto ret = FieldVector<double, d>({1. + rand1, 0. + rand2});
                                                           ret *= (phi_initial.evaluate(x, param) + 1.) / 2.;
                                                           return ret;
                                                         },
                                                         /*name=*/"P_initial");

    interpolate(phi_initial, phi);
    interpolate(P_initial, P);

    // implicit Euler timestepping
    double t = 0;
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    double next_save_time = t + write_step > t_end ? t_end : t + write_step;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    write_files(true, u, p, P, Pnat, phi, phinat, mu, filename, 0, t);

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      double max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      solve_navier_stokes(u, p, P, Pnat, phi, phinat, mu, Re, Fa, xi, dirichlet_boundary_info);
      solve_ofield(u, p, P, Pnat, phi, phinat, mu, xi, kappa, c_1, Pa, beta);
      solve_pfield(u, p, P, Pnat, phi, phinat, mu, gamma, c_1, Pa, beta, epsilon);

      t += actual_dt;

      // check if data should be written in this timestep (and write)
      if (Dune::XT::Common::FloatCmp::ge(t, next_save_time)) {
        write_files(true, u, p, P, Pnat, phi, phinat, mu, filename, save_step_counter, t);
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
