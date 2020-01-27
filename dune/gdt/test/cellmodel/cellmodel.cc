// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <chrono>
#include <cstdlib>

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-market.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/solver/istl/preconditioners.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/symmetrized-laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/gradient-value.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/norms.hh>

#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>

#include <Eigen/IterativeLinearSolvers>

#include "fgmres.hh"
#include "linearsolvers.hh"
#include "cellmodel.hh"

using namespace Dune;
using namespace Dune::GDT;

template <class DiscreteFunc>
static std::string get_filename(const std::string& prefix, const DiscreteFunc& func, const size_t step)
{
  return prefix + "_" + func.name() + "_" + XT::Common::to_string(step) + ".txt";
}

template <class DiscreteFunc>
static std::string
get_rankfilename(const std::string& prefix, const DiscreteFunc& func, const size_t step, const int rank)
{
  return prefix + "_rank_" + XT::Common::to_string(rank) + "_" + func.name() + "_" + XT::Common::to_string(step)
         + ".txt";
}

template <class DiscreteFunc>
static void write_to_textfile(const DiscreteFunc& func, const std::string& prefix, const size_t step, const double t)
{
  static const size_t d = DiscreteFunc::d;
  const auto& grid_view = func.space().grid_view();
  const auto local_func = func.local_function();
  // write one file per MPI rank
  std::ofstream rankfile(get_rankfilename(prefix, func, step, grid_view.comm().rank()));
  for (const auto& entity : elements(grid_view, Partitions::interiorBorder)) {
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

CellModelSolver::CellModelSolver(const std::string testcase,
                                 const double t_end,
                                 const unsigned int num_elements_x,
                                 const unsigned int num_elements_y,
                                 const bool use_tbb,
                                 const double Be, // bending capillary number, ratio of viscous forces to bending forces
                                 const double Ca, // capillary number, ratio of viscous forces to surface tension forces
                                 const double Pa, // polarization elasticity number
                                 const double Re, // Reynolds number
                                 const double Fa, // active force number
                                 const double xi, // alignment of P with the flow, > 0 for rod-like cells and < 0 for
                                                  // oblate ones
                                 const double kappa, // eta_rot/eta, scaling factor between rotational and dynamic
                                                     // viscosity
                                 const double c_1, // double well shape parameter
                                 const double beta, // alignment of P with the boundary of cell
                                 const double gamma, // phase field mobility coefficient
                                 const double epsilon, // phase field parameter
                                 const double In, // interaction parameter
                                 const CellModelLinearSolverType pfield_solver_type,
                                 const CellModelMassMatrixSolverType pfield_mass_matrix_solver_type,
                                 const CellModelLinearSolverType ofield_solver_type,
                                 const CellModelMassMatrixSolverType ofield_mass_matrix_solver_type,
                                 const double outer_reduction,
                                 const int outer_restart,
                                 const int outer_verbose,
                                 const double inner_reduction,
                                 const int inner_maxit,
                                 const int inner_verbose,
                                 const bool linearize,
                                 const double pol_order)
  : lower_left_(get_lower_left(testcase))
  , upper_right_(get_upper_right(testcase))
  , t_end_(t_end)
  , t_(0.)
  , use_tbb_(use_tbb)
  , Re_(Re)
  , Fa_inv_(1. / Fa)
  , xi_(xi)
  , kappa_(kappa)
  , c_1_(c_1)
  , Pa_(Pa)
  , last_pfield_Pa_(Pa_)
  , last_ofield_Pa_(Pa_)
  , beta_(beta)
  , gamma_(gamma)
  , Be_(Be)
  , Ca_(Ca)
  , epsilon_(epsilon)
  , In_(In)
  , vol_domain_((upper_right_[0] - lower_left_[0]) * (upper_right_[1] - lower_left_[1]))
  , num_cells_(get_num_cells(testcase))
  , linearize_(linearize)
  // do a global refine once, this makes simplicial grids look more symmetric
  , grid_(XT::Grid::make_cube_grid<G>(lower_left_, upper_right_, {num_elements_x, num_elements_y}, 1))
  , nonperiodic_grid_view_(grid_.leaf_view())
  , grid_view_(nonperiodic_grid_view_, std::bitset<d>(get_periodic_directions(testcase)))
  , u_space_(make_continuous_lagrange_space<d>(grid_view_, pol_order))
  , p_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order - 1))
  , phi_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order))
  , size_u_(u_space_.mapper().size())
  , size_p_(p_space_.mapper().size())
  , size_phi_(phi_space_.mapper().size())
  , num_mutexes_u_(use_tbb ? 0 : size_u_ / 100)
  , num_mutexes_ofield_(use_tbb ? 0 : size_u_ / 100)
  , num_mutexes_pfield_(use_tbb ? 0 : size_phi_ / 100)
  , stokes_vector_(size_u_ + size_p_, 0., 0)
  , ofield_vectors_(num_cells_, VectorType(2 * size_u_, 0., 0))
  , pfield_vectors_(num_cells_, VectorType(3 * size_phi_, 0., 0))
  , u_view_(stokes_vector_, 0, size_u_)
  , p_view_(stokes_vector_, size_u_, size_u_ + size_p_)
  , u_(u_space_, u_view_, "u")
  , p_(p_space_, p_view_, "p")
  , S_stokes_(size_u_ + size_p_, size_u_ + size_p_, create_stokes_pattern(u_space_, p_space_), 100)
  , A_stokes_(S_stokes_, 0, size_u_, 0, size_u_)
  , B_stokes_(S_stokes_, 0, size_u_, size_u_, size_u_ + size_p_)
  , BT_stokes_(S_stokes_, size_u_, size_u_ + size_p_, 0, size_u_)
  , M_p_stokes_(size_p_, size_p_, make_element_sparsity_pattern(p_space_, p_space_, grid_view_), 100)
  , A_stokes_op_(std::make_shared<MatrixOperator<MatrixViewType, PGV, d>>(grid_view_, u_space_, u_space_, A_stokes_))
  , stokes_rhs_vector_(size_u_ + size_p_, 0., num_mutexes_u_)
  , stokes_f_vector_(stokes_rhs_vector_, 0, size_u_)
  , stokes_g_vector_(stokes_rhs_vector_, size_u_, size_u_ + size_p_)
  , p_basis_integrated_vector_(size_p_)
  , u_dirichlet_constraints_(make_dirichlet_constraints(u_space_, boundary_info_))
  , phi_dirichlet_constraints_(make_dirichlet_constraints(phi_space_, boundary_info_))
  , ofield_submatrix_pattern_(make_element_sparsity_pattern(u_space_, u_space_, grid_view_))
  , M_ofield_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , A_ofield_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_elliptic_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_linear_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_nonlinear_part_(size_u_, size_u_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , S_schur_ofield_linear_part_(size_u_, size_u_, ofield_submatrix_pattern_, 0)
  , M_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, M_ofield_))
  , A_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, A_ofield_))
  , C_ofield_linear_part_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, C_ofield_linear_part_))
  , C_ofield_nonlinear_part_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, C_ofield_nonlinear_part_))
  , ofield_jac_linear_op_(M_ofield_, A_ofield_, C_ofield_linear_part_, kappa_, this)
  , ofield_solver_(kappa_,
                   M_ofield_,
                   A_ofield_,
                   C_ofield_linear_part_,
                   C_ofield_nonlinear_part_,
                   S_schur_ofield_linear_part_,
                   ofield_solver_type,
                   ofield_mass_matrix_solver_type,
                   ofield_submatrix_pattern_,
                   num_cells_,
                   outer_reduction,
                   outer_restart,
                   outer_verbose,
                   inner_reduction,
                   inner_maxit,
                   inner_verbose)
  , ofield_rhs_vector_(2 * size_u_, 0., num_mutexes_ofield_)
  , ofield_f_vector_(ofield_rhs_vector_, 0, size_u_)
  , ofield_g_vector_(ofield_rhs_vector_, size_u_, 2 * size_u_)
  , stokes_solver_(std::make_shared<LUSolverType>())
  , ofield_tmp_vec_(2 * size_u_, 0., 0)
  , ofield_tmp_vec2_(2 * size_u_, 0., 0)
  , ofield_deim_input_dofs_(num_cells_)
  , Pnat_deim_input_dofs_begin_(num_cells_)
  , ofield_deim_output_dofs_(num_cells_)
  , ofield_deim_unique_output_dofs_(num_cells_)
  , P_deim_output_dofs_(num_cells_)
  , Pnat_deim_output_dofs_(num_cells_)
  , ofield_deim_entities_(num_cells_)
  , pfield_submatrix_pattern_(make_element_sparsity_pattern(phi_space_, phi_space_, grid_view_))
  , M_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , D_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , M_ell_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , M_nonlin_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , G_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , A_boundary_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , D_pfield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, D_pfield_))
  , G_pfield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, G_pfield_))
  , M_nonlin_pfield_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, M_nonlin_pfield_))
  , pfield_jac_linear_op_(M_pfield_,
                          D_pfield_,
                          M_ell_pfield_,
                          A_boundary_pfield_,
                          phi_dirichlet_constraints_,
                          gamma_,
                          epsilon_,
                          Be_,
                          this)
  , pfield_solver_(gamma_,
                   epsilon_,
                   Be_,
                   Ca_,
                   M_pfield_,
                   M_ell_pfield_,
                   D_pfield_,
                   G_pfield_,
                   M_nonlin_pfield_,
                   A_boundary_pfield_,
                   pfield_solver_type,
                   pfield_mass_matrix_solver_type,
                   phi_dirichlet_constraints_.dirichlet_DoFs(),
                   pfield_submatrix_pattern_,
                   num_cells_,
                   outer_reduction,
                   outer_restart,
                   outer_verbose,
                   inner_reduction,
                   inner_maxit,
                   inner_verbose)
  , pfield_rhs_vector_(3 * size_phi_, 0., num_mutexes_pfield_)
  , pfield_g_vector_(pfield_rhs_vector_, 0, size_phi_)
  , pfield_h_vector_(pfield_rhs_vector_, size_phi_, 2 * size_phi_)
  , pfield_f_vector_(pfield_rhs_vector_, 2 * size_phi_, 3 * size_phi_)
  , pfield_deim_input_dofs_(num_cells_)
  , phinat_deim_input_dofs_begin_(num_cells_)
  , mu_deim_input_dofs_begin_(num_cells_)
  , pfield_deim_output_dofs_(num_cells_)
  , pfield_deim_unique_output_dofs_(num_cells_)
  , phi_deim_output_dofs_(num_cells_)
  , phinat_deim_output_dofs_(num_cells_)
  , mu_deim_output_dofs_(num_cells_)
  , both_mu_and_phi_deim_output_dofs_(num_cells_)
  , pfield_deim_entities_(num_cells_)
  , pfield_tmp_vec_(3 * size_phi_, 0., 0)
  , pfield_tmp_vec2_(3 * size_phi_, 0., 0)
  , phi_tmp_vec_(size_phi_, 0., 0)
  , u_tmp_vec_(size_u_, 0., 0)
  , u_tmp_(u_space_)
  , u_tmp_local_(std::make_shared<PerThreadVectorLocalFunc>())
  , P_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
  , Pnat_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
  , phi_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
  , phinat_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
  , mu_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
{

  /************************** create and project initial values*****************************************
   ************************** we only need initial values for P and phi ********************************
   ************************** mu_initial is only needed if linearization is used ***********************/

  std::shared_ptr<const XT::Functions::FunctionInterface<d, d>> u_initial_func;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> phi_initial_funcs;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> mu_initial_funcs;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d, d>>> P_initial_funcs;

  if (testcase == "single_cell") {
    // Initially, cell is circular with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 160] \times [0, 40].
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(80, 20) - x|.
    FieldVector<double, d> center{upper_right_[0] / 2., upper_right_[1] / 2.};
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
    FieldVector<double, d> center1{15, 15};
    FieldVector<double, d> center2{35, 35};
    auto r1 = [center1](const auto& xr) { return 4.0 - (center1 - xr).two_norm(); };
    auto r2 = [center2](const auto& xr) { return 4.0 - (center2 - xr).two_norm(); };
    const XT::Functions::GenericFunction<d> phi1_initial(
        50,
        /*evaluate=*/
        [r = r1, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2.) * epsilon)); },
        /*name=*/"phi1_initial");
    const XT::Functions::GenericFunction<d> mu1_initial(50,
                                                        /*evaluate=*/
                                                        [phi1_initial, epsilon](const auto& x, const auto& param) {
                                                          // TODO: add approximation of laplacian term
                                                          const auto phi = phi1_initial.evaluate(x, param);
                                                          return 1. / epsilon * (std::pow(phi, 3) - phi);
                                                        },
                                                        /*name=*/"mu1_initial");
    const XT::Functions::GenericFunction<d> phi2_initial(
        50,
        /*evaluate=*/
        [r = r2, epsilon](const auto& x, const auto& /*param*/) { return std::tanh(r(x) / (std::sqrt(2.) * epsilon)); },
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

  /*************************************************************************************************
   ******************************* create variables, set initial values ****************************
   *************************************************************************************************/

  // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
  // conditions for the other variables
  const XT::Functions::ConstantFunction<d> minus_one(-1.);
  XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
  default_interpolation(*u_initial_func, u_);
  // create system and temporary vectors, DiscreteFunctions, etc.
  for (size_t kk = 0; kk < num_cells_; kk++) {
    P_view_.emplace_back(ofield_vectors_[kk], 0, size_u_);
    Pnat_view_.emplace_back(ofield_vectors_[kk], size_u_, 2 * size_u_);
    phi_view_.emplace_back(pfield_vectors_[kk], 0, size_phi_);
    phinat_view_.emplace_back(pfield_vectors_[kk], size_phi_, 2 * size_phi_);
    mu_view_.emplace_back(pfield_vectors_[kk], 2 * size_phi_, 3 * size_phi_);
    const auto kk_str = XT::Common::to_string(kk);
    P_.emplace_back(make_discrete_function(u_space_, P_view_[kk], "P_" + kk_str));
    Pnat_.emplace_back(make_discrete_function(u_space_, Pnat_view_[kk], "Pnat_" + kk_str));
    phi_.emplace_back(make_discrete_function(phi_space_, phi_view_[kk], "phi_" + kk_str));
    phinat_.emplace_back(make_discrete_function(phi_space_, phinat_view_[kk], "phinat_" + kk_str));
    mu_.emplace_back(make_discrete_function(phi_space_, mu_view_[kk], "mu_" + kk_str));
    P_tmp_.emplace_back(u_space_);
    Pnat_tmp_.emplace_back(u_space_);
    phi_tmp_.emplace_back(phi_space_);
    phinat_tmp_.emplace_back(phi_space_);
    mu_tmp_.emplace_back(phi_space_);
    default_interpolation(*phi_initial_funcs[kk], phi_[kk]);
    boundary_interpolation(minus_one, phi_[kk], all_dirichlet_boundary_info, XT::Grid::DirichletBoundary{});
    default_interpolation(*mu_initial_funcs[kk], mu_[kk]);
    default_interpolation(*P_initial_funcs[kk], P_[kk]);
  }

  /*************************************************************************************************
   *************************************** Stokes **************************************************
   *************************************************************************************************/
  if (Re_ > 1e-2)
    DUNE_THROW(NotImplemented, "No Navier-Stokes solver implemented yet!");

  MatrixOperator<MatrixViewType, PGV, 1, 1, d> B_stokes_op(grid_view_, p_space_, u_space_, B_stokes_);
  MatrixOperator<MatrixType, PGV, 1, 1, 1> M_p_stokes_op(grid_view_, p_space_, p_space_, M_p_stokes_);
  // calculate A_{ij} as \int \nabla v_i \nabla v_j
  // A_stokes_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalSymmetricEllipticIntegrand<E>(1.)));
  A_stokes_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>()));
  // calculate B_{ij} as \int \nabla p_i div(v_j)
  B_stokes_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
      LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
  M_p_stokes_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>()));

  auto p_basis_integrated_functional = make_vector_functional(p_space_, p_basis_integrated_vector_);
  const XT::Functions::ConstantGridFunction<E> one_function(1);
  p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), one_function)));
  B_stokes_op.append(p_basis_integrated_functional);

  // Dirichlet constrainst for u
  A_stokes_op_->append(u_dirichlet_constraints_);
  // assemble everything
  A_stokes_op_->append(B_stokes_op);
  A_stokes_op_->append(M_p_stokes_op);
  A_stokes_op_->assemble(use_tbb_);

  // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the p_size_-th row of
  // [A B; B^T 0] to the unit vector.
  const size_t dof_index = 0;
  S_stokes_.set_entry(size_u_ + dof_index, size_u_ + dof_index, 1.);
  B_stokes_.clear_col(dof_index);
  stokes_g_vector_.set_entry(dof_index, 0.);

  u_dirichlet_constraints_.apply(A_stokes_, false, true);
  for (const auto& DoF : u_dirichlet_constraints_.dirichlet_DoFs())
    B_stokes_.clear_row(DoF);

  // Set B^T
  const auto B_pattern = B_stokes_.pattern();
  for (size_t ii = 0; ii < size_u_; ii++)
    for (const auto& jj : B_pattern.inner(ii))
      BT_stokes_.set_entry(jj, ii, B_stokes_.get_entry(ii, jj));

  S_colmajor_ = S_stokes_.backend();
  stokes_solver_->analyzePattern(S_colmajor_);
  stokes_solver_->factorize(S_colmajor_);

  /*************************************************************************************************
   ************************************ Orientationfield *******************************************
   *************************************************************************************************/
  // calculate M_{ij} as \int \psi_i phi_j
  M_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(1.)));
  C_ofield_elliptic_part_ *= 0.;
  MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, u_space_, u_space_, C_ofield_elliptic_part_);
  ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
  M_ofield_op_->append(ofield_elliptic_op);
  M_ofield_op_->assemble(use_tbb_);
  ofield_solver_.setup();

  /*************************************************************************************************
   **************************************** Phasefield *********************************************
   *************************************************************************************************/

  MatrixOperator<MatrixType, PGV, 1> M_pfield_op(grid_view_, phi_space_, phi_space_, M_pfield_);
  M_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(1.)));
  MatrixOperator<MatrixType, PGV, 1> M_ell_pfield_op(grid_view_, phi_space_, phi_space_, M_ell_pfield_);
  M_ell_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
  MatrixOperator<MatrixType, PGV, 1> A_boundary_pfield_op(grid_view_, phi_space_, phi_space_, A_boundary_pfield_);
  A_boundary_pfield_op.append(
      LocalIntersectionIntegralBilinearForm<PI, 1>(LocalBoundaryIntersectionGradientValueIntegrand<PI>(1.)),
      {},
      XT::Grid::ApplyOn::NonPeriodicBoundaryIntersections<PGV>());
  M_pfield_op.append(phi_dirichlet_constraints_);
  M_pfield_op.append(M_ell_pfield_op);
  M_pfield_op.append(A_boundary_pfield_op);
  M_pfield_op.assemble(use_tbb_);
  pfield_solver_.setup();
} // constructor

size_t CellModelSolver::num_cells() const
{
  return num_cells_;
}

bool CellModelSolver::linear() const
{
  return linearize_;
}

bool CellModelSolver::finished() const
{
  return XT::Common::FloatCmp::eq(t_end_, t_);
}

//******************************************************************************************************************
//********************************* Solve methods for whole system of equations ************************************
//******************************************************************************************************************

// Solves whole system of equations using the values stored in stokes_vector_, ofield_vectors_ and pfield_vectors_  as
// initial values Returns the whole trajectory, i.e., ret[i] contains the results in the i-th timestep. The first
// num_cells_ entries of ret[i] correspond to the phasefield for each cell, the next num_cells entries are the
// orientation field vectors, and the last one is the stokes vector. dt: Time step length. write: Whether to write
// .vtu and .txt files. write_step: Time interval at which results should be written. If negative, all steps are
// written. Ignored if write = false. filename: Prefix for .vtu and .txt files. Ignored if write = false. subsampling:
// Whether to use subsampling for visualization. Ignored if write = false.
std::vector<std::vector<CellModelSolver::VectorType>> CellModelSolver::solve(
    const double dt, const bool write, const double write_step, const std::string filename, const bool subsampling)
{
  std::vector<std::vector<VectorType>> ret(1 + 2 * num_cells_);
  for (size_t kk = 0; kk < num_cells_; ++kk) {
    ret[kk].push_back(pfield_vectors_[kk]);
    ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
  }
  ret[2 * num_cells_].push_back(stokes_vector_);
  // implicit Euler timestepping
  if (write)
    visualize(filename, 0, t_, subsampling);
  assert(Dune::XT::Common::FloatCmp::ge(t_end_, t_));
  double next_save_time = t_ + write_step > t_end_ ? t_end_ : t_ + write_step;
  size_t save_step_counter = 1;

  while (Dune::XT::Common::FloatCmp::lt(t_, t_end_)) {
    double max_dt = dt;
    // match saving times and t_end_ exactly
    if (Dune::XT::Common::FloatCmp::gt(t_ + dt, t_end_))
      max_dt = t_end_ - t_;
    double actual_dt = std::min(dt, max_dt);

    // do a timestep
    std::cout << "Current time: " << t_ << std::endl;
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      prepare_pfield_operator(dt, kk);
      ret[kk].push_back(apply_inverse_pfield_operator(ret[kk].back(), kk));
      set_pfield_vec(kk, ret[kk].back());
      std::cout << "Pfield " << kk << " done" << std::endl;
      prepare_ofield_operator(dt, kk);
      ret[num_cells_ + kk].push_back(apply_inverse_ofield_operator(ret[num_cells_ + kk].back(), kk));
      set_ofield_vec(kk, ret[num_cells_ + kk].back());
      std::cout << "Ofield " << kk << " done" << std::endl;
    }

    // stokes system
    prepare_stokes_operator();
    ret[2 * num_cells_].push_back(apply_inverse_stokes_operator());
    set_stokes_vec(ret[2 * num_cells_].back());
    std::cout << "Stokes done" << std::endl;

    t_ += actual_dt;

    // check if data should be written in this timestep (and write)
    if (write) {
      if (write_step < 0. || Dune::XT::Common::FloatCmp::ge(t_, next_save_time)) {
        visualize(filename, save_step_counter, t_, subsampling);
        next_save_time += write_step;
        ++save_step_counter;
      }
    }
  } // while (t_ < t_end_)
  return ret;
}

// Like solve, but only computes and returns the next n timesteps
std::vector<std::vector<CellModelSolver::VectorType>> CellModelSolver::next_n_timesteps(const size_t n, const double dt)
{
  std::vector<std::vector<VectorType>> ret(2 * num_cells_ + 1);
  size_t count = 0;
  if (XT::Common::is_zero(t_)) {
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      ret[kk].push_back(pfield_vectors_[kk]);
      ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
    }
    ret[2 * num_cells_].push_back(stokes_vector_);
    // Hack to avoid adding initial_values twice
    ++count;
    t_ = 1e-100;
  }
  // Undo hack to avoid adding initial_values twice
  if (XT::Common::is_zero(t_ - 1e-100) && count == 0)
    t_ = 0.;

  assert(Dune::XT::Common::FloatCmp::ge(t_end_, t_));

  // implicit Euler timestepping
  while (Dune::XT::Common::FloatCmp::lt(t_, t_end_) && count < n) {
    double max_dt = dt;
    // match saving times and t_end_ exactly
    if (Dune::XT::Common::FloatCmp::gt(t_ + dt, t_end_))
      max_dt = t_end_ - t_;
    double actual_dt = std::min(dt, max_dt);

    // do a timestep
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      prepare_pfield_operator(dt, kk);
      pfield_vectors_[kk] = apply_inverse_pfield_operator(pfield_vectors_[kk], kk);
      ret[kk].push_back(pfield_vectors_[kk]);
      // std::cout << "Pfield " << kk << " done" << std::endl;
      prepare_ofield_operator(dt, kk);
      ofield_vectors_[kk] = apply_inverse_ofield_operator(ofield_vectors_[kk], kk);
      ret[num_cells_ + kk].push_back(ofield_vectors_[kk]);
      // std::cout << "Ofield " << kk << " done" << std::endl;
    }

    // stokes system
    prepare_stokes_operator();
    stokes_vector_ = apply_inverse_stokes_operator();
    ret[2 * num_cells_].push_back(stokes_vector_);
    // std::cout << "Stokes done" << std::endl;

    ++count;
    t_ += actual_dt;
  } // while (t_ < t_end_)
  return ret;
}

//******************************************************************************************************************
//********************************* Product operators (mass matrix application) ************************************
//******************************************************************************************************************

// applies the pfield mass matrix to phi, phinat, mu
// To calculate the sum of the squared L2 products of phi, phinat and mu, calculate the inner product of the result
// with vec.
CellModelSolver::VectorType CellModelSolver::apply_pfield_product_operator(const VectorType& vec) const
{
  VectorType ret(3 * size_phi_);
  ConstVectorViewType phi_view(vec, 0, size_phi_);
  ConstVectorViewType phinat_view(vec, size_phi_, 2 * size_phi_);
  ConstVectorViewType mu_view(vec, 2 * size_phi_, 3 * size_phi_);
  VectorViewType phi_ret_view(ret, 0, size_phi_);
  VectorViewType phinat_ret_view(ret, size_phi_, 2 * size_phi_);
  VectorViewType mu_ret_view(ret, 2 * size_phi_, 3 * size_phi_);
  M_pfield_.mv(phi_view, phi_ret_view);
  M_pfield_.mv(phinat_view, phinat_ret_view);
  M_pfield_.mv(mu_view, mu_ret_view);
  return ret;
}

// applies the ofield mass matrix to P, Pnat
// To calculate the sum of the squared L2 products of P and Pnat, calculate the inner product of the result with vec.
CellModelSolver::VectorType CellModelSolver::apply_ofield_product_operator(const VectorType& vec) const
{
  VectorType ret(2 * size_u_);
  ConstVectorViewType P_view(vec, 0, size_u_);
  ConstVectorViewType Pnat_view(vec, size_u_, 2 * size_u_);
  VectorViewType P_ret_view(ret, 0, size_u_);
  VectorViewType Pnat_ret_view(ret, size_u_, 2 * size_u_);
  M_ofield_.mv(P_view, P_ret_view);
  M_ofield_.mv(Pnat_view, Pnat_ret_view);
  return ret;
}

// applies the ofield mass matrix to P, Pnat
// To calculate the sum of the squared L2 products of P and Pnat, calculate the inner product of the result with vec.
CellModelSolver::VectorType CellModelSolver::apply_stokes_product_operator(const VectorType& vec) const
{
  VectorType ret(size_u_ + size_p_);
  ConstVectorViewType u_view(vec, 0, size_u_);
  ConstVectorViewType p_view(vec, size_u_, size_u_ + size_p_);
  VectorViewType u_ret_view(ret, 0, size_u_);
  VectorViewType p_ret_view(ret, size_u_, size_u_ + size_p_);
  // The Orientation field variables and u have the same basis so use M_ofield_
  M_ofield_.mv(u_view, u_ret_view);
  M_p_stokes_.mv(p_view, p_ret_view);
  return ret;
}

//******************************************************************************************************************
//*****************************************  Visualization   *******************************************************
//******************************************************************************************************************

// Visualizes given vector as phasefield finite element vector
void CellModelSolver::visualize_pfield(const std::string& filename, const VectorType& vec, const bool subsampling) const
{
  auto vtk_writer = phi_[0].create_vtkwriter(phi_space_.grid_view(), subsampling);
  const ConstVectorViewType phi_vec(vec, 0, size_phi_);
  const ConstVectorViewType phinat_vec(vec, size_phi_, 2 * size_phi_);
  const ConstVectorViewType mu_vec(vec, 2 * size_phi_, 3 * size_phi_);
  const auto phi_func = make_discrete_function(phi_space_, phi_vec, "phi");
  const auto phinat_func = make_discrete_function(phi_space_, phinat_vec, "phinat");
  const auto mu_func = make_discrete_function(phi_space_, mu_vec, "mu");
  phi_func.add_to_vtkwriter(*vtk_writer);
  phinat_func.add_to_vtkwriter(*vtk_writer);
  mu_func.add_to_vtkwriter(*vtk_writer);
  phi_[0].write_visualization(*vtk_writer, filename);
} // void visualize_pfield(...)

// Visualizes given vector as orientation field finite element vector
void CellModelSolver::visualize_ofield(const std::string& filename, const VectorType& vec, const bool subsampling) const
{
  auto vtk_writer = P_[0].create_vtkwriter(u_space_.grid_view(), subsampling);
  const ConstVectorViewType P_vec(vec, 0, size_u_);
  const ConstVectorViewType Pnat_vec(vec, size_u_, 2 * size_u_);
  const auto P_func = make_discrete_function(u_space_, P_vec, "P");
  const auto Pnat_func = make_discrete_function(u_space_, Pnat_vec, "Pnat");
  P_func.add_to_vtkwriter(*vtk_writer);
  Pnat_func.add_to_vtkwriter(*vtk_writer);
  P_[0].write_visualization(*vtk_writer, filename);
} // void visualize_ofield(...)

// Visualizes given vector as stokes finite element vector
void CellModelSolver::visualize_stokes(const std::string& filename, const VectorType& vec, const bool subsampling) const
{
  auto vtk_writer = u_.create_vtkwriter(u_.space().grid_view(), subsampling);
  const ConstVectorViewType u_vec(vec, 0, size_u_);
  const ConstVectorViewType p_vec(vec, size_u_, size_u_ + size_p_);
  const auto u_func = make_discrete_function(u_space_, u_vec, "u");
  const auto p_func = make_discrete_function(p_space_, p_vec, "p");
  u_.add_to_vtkwriter(*vtk_writer);
  p_func.add_to_vtkwriter(*vtk_writer);
  u_.write_visualization(*vtk_writer, filename);
} // void visualize_stokes(...)

// Visualizes variables currently stored in this class.
// If txt = true, also writes textfiles containing the values.
void CellModelSolver::visualize(const std::string& prefix,
                                const size_t step,
                                const double t,
                                const bool subsampling,
                                const bool vtu,
                                const bool txt) const
{
  auto vtk_writer = u_.create_vtkwriter(u_.space().grid_view(), subsampling);
  std::string postfix = "_" + XT::Common::to_string(step);
  if (vtu) {
    u_.add_to_vtkwriter(*vtk_writer);
    p_.add_to_vtkwriter(*vtk_writer);
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      // std::cout << "phi l2_norm: " << l2_norm(phi_[kk].space().grid_view(), phi_[kk]) << std::endl;
      // std::cout << "phinat l2_norm: " << l2_norm(phinat_[kk].space().grid_view(), phinat_[kk]) << std::endl;
      // std::cout << "mu l2_norm: " << l2_norm(mu_[kk].space().grid_view(), mu_[kk]) << std::endl;
      P_[kk].add_to_vtkwriter(*vtk_writer);
      Pnat_[kk].add_to_vtkwriter(*vtk_writer);
      phi_[kk].add_to_vtkwriter(*vtk_writer);
      phinat_[kk].add_to_vtkwriter(*vtk_writer);
      mu_[kk].add_to_vtkwriter(*vtk_writer);
      phi_[kk].add_gradient_to_vtkwriter(*vtk_writer);
      phinat_[kk].add_gradient_to_vtkwriter(*vtk_writer);
      mu_[kk].add_gradient_to_vtkwriter(*vtk_writer);
    }
    u_.write_visualization(*vtk_writer, prefix + postfix);
  } // if (vtu)
  if (txt) {
    write_to_textfile(u_, prefix, step, t);
    write_to_textfile(p_, prefix, step, t);
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      write_to_textfile(P_[kk], prefix, step, t);
      write_to_textfile(Pnat_[kk], prefix, step, t);
      write_to_textfile(phi_[kk], prefix, step, t);
      write_to_textfile(phinat_[kk], prefix, step, t);
      write_to_textfile(mu_[kk], prefix, step, t);
    } // kk
  } // if (txt)
} // void visualize(...)

//******************************************************************************************************************
//*******************************  Methods to get and set variable values   ****************************************
//******************************************************************************************************************

// Sets stokes vector to stokes_vec
void CellModelSolver::set_stokes_vec(const VectorType& stokes_vec)
{
  DUNE_THROW_IF(
      stokes_vec.size() != size_u_ + size_p_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
  stokes_vector_ = stokes_vec;
}

// Sets orientation field vector belonging to cell to pfield_vec
void CellModelSolver::set_ofield_vec(const size_t cell, const VectorType& ofield_vec)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(ofield_vec.size() != 2 * size_u_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
  ofield_vectors_[cell] = ofield_vec;
}

// Sets phasefield vector belonging to cell to pfield_vec
void CellModelSolver::set_pfield_vec(const size_t cell, const VectorType& pfield_vec)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(pfield_vec.size() != 3 * size_phi_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
  pfield_vectors_[cell] = pfield_vec;
}

// Get stokes finite element vector
const CellModelSolver::VectorType& CellModelSolver::stokes_vec()
{
  return stokes_vector_;
}

// Get orientation field finite element vector belonging to cell
const CellModelSolver::VectorType& CellModelSolver::ofield_vec(const size_t cell)
{
  return ofield_vectors_[cell];
}

// Get phase field finite element vector belonging to cell
const CellModelSolver::VectorType& CellModelSolver::pfield_vec(const size_t cell)
{
  return pfield_vectors_[cell];
}

//******************************************************************************************************************
//****** Prepare methods (calculate everything that is linear for the respective operator, but depends on **********
//****** the values of other variables, so cannot be computed once and for all in the constructor )       **********
//******************************************************************************************************************

void CellModelSolver::prepare_stokes_operator()
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  for (size_t kk = 0; kk < num_cells_; kk++) {
    phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
    phinat_tmp_[kk].dofs().vector() = phinat_[kk].dofs().vector();
    P_tmp_[kk].dofs().vector() = P_[kk].dofs().vector();
    Pnat_tmp_[kk].dofs().vector() = Pnat_[kk].dofs().vector();
  }
  assemble_stokes_rhs();
}

void CellModelSolver::prepare_ofield_operator(const double dt, const size_t cell, const bool restricted)
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
  phi_tmp_[cell].dofs().vector() = phi_[cell].dofs().vector();
  dt_ = dt;
  ofield_jac_linear_op_.prepare(dt, cell, restricted);
  assemble_ofield_rhs(dt, cell);
  assemble_ofield_linear_jacobian(dt, cell);
}

void CellModelSolver::prepare_pfield_operator(const double dt, const size_t cell, const bool restricted)
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
  for (size_t kk = 0; kk < num_cells_; kk++) {
    phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
    mu_tmp_[kk].dofs().vector() = mu_[kk].dofs().vector();
  }
  assemble_pfield_rhs(dt, cell, restricted);
  assemble_pfield_linear_jacobian(dt, cell, restricted);
  pfield_jac_linear_op_.prepare(dt, cell, restricted);
  dt_ = dt;
}

void CellModelSolver::compute_restricted_ofield_dofs(const std::vector<size_t>& output_dofs, const size_t cell)
{
  if (!ofield_deim_output_dofs_[cell] || *ofield_deim_output_dofs_[cell] != output_dofs) {
    const auto& pattern = ofield_solver_.system_matrix_pattern(ofield_submatrix_pattern_);
    // We need to keep the original output_dofs which is unordered and may contain duplicates, as the restricted
    // operator will return exactly these dofs. For computations, however, we often need unique dofs.
    ofield_deim_output_dofs_[cell] = std::make_shared<std::vector<size_t>>(output_dofs);
    auto& unique_output_dofs = ofield_deim_unique_output_dofs_[cell];
    unique_output_dofs = output_dofs;
    std::sort(unique_output_dofs.begin(), unique_output_dofs.end());
    unique_output_dofs.erase(std::unique(unique_output_dofs.begin(), unique_output_dofs.end()),
                             unique_output_dofs.end());
    // sort output into dofs belonging to phi, phinat and mu
    auto& P_output_dofs = P_deim_output_dofs_[cell];
    auto& Pnat_output_dofs = Pnat_deim_output_dofs_[cell];
    P_output_dofs.clear();
    Pnat_output_dofs.clear();
    for (const auto& dof : unique_output_dofs) {
      if (dof < size_u_)
        P_output_dofs.push_back(dof);
      else
        Pnat_output_dofs.push_back(dof);
    }
    for (auto& dof : Pnat_output_dofs)
      dof -= size_u_;
    // get input dofs corresponding to output dofs
    auto& input_dofs = ofield_deim_input_dofs_[cell];
    input_dofs.clear();
    for (const auto& dof : unique_output_dofs) {
      const auto& new_input_dofs = pattern.inner(dof);
      input_dofs.insert(input_dofs.end(), new_input_dofs.begin(), new_input_dofs.end());
    }
    // sort and remove duplicate entries
    std::sort(input_dofs.begin(), input_dofs.end());
    input_dofs.erase(std::unique(input_dofs.begin(), input_dofs.end()), input_dofs.end());
    Pnat_deim_input_dofs_begin_[cell] =
        std::lower_bound(input_dofs.begin(), input_dofs.end(), size_u_) - input_dofs.begin();
    // store all entities that contain an output dof
    const auto& mapper = u_space_.mapper();
    DynamicVector<size_t> global_indices;
    ofield_deim_entities_[cell].clear();
    for (const auto& entity : Dune::elements(grid_view_)) {
      mapper.global_indices(entity, global_indices);
      maybe_add_entity(entity, global_indices, *ofield_deim_output_dofs_[cell], ofield_deim_entities_[cell], size_u_);
    } // entities
  } // if (not already computed)
} // void compute_restricted_pfield_dofs(...)

void CellModelSolver::compute_restricted_pfield_dofs(const std::vector<size_t>& output_dofs, const size_t cell)
{
  if (!pfield_deim_output_dofs_[cell] || *pfield_deim_output_dofs_[cell] != output_dofs) {
    const auto& pattern = pfield_solver_.system_matrix_pattern(pfield_submatrix_pattern_);
    // We need to keep the original output_dofs which is unordered and may contain duplicates, as the restricted
    // operator will return exactly these dofs. For computations, however, we often need unique dofs.
    pfield_deim_output_dofs_[cell] = std::make_shared<std::vector<size_t>>(output_dofs);
    auto& unique_output_dofs = pfield_deim_unique_output_dofs_[cell];
    unique_output_dofs = output_dofs;
    std::sort(unique_output_dofs.begin(), unique_output_dofs.end());
    unique_output_dofs.erase(std::unique(unique_output_dofs.begin(), unique_output_dofs.end()),
                             unique_output_dofs.end());
    // sort output into dofs belonging to phi, phinat and mu
    auto& phi_output_dofs = phi_deim_output_dofs_[cell];
    auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
    auto& mu_output_dofs = mu_deim_output_dofs_[cell];
    auto& phinat_mu_output_dofs = both_mu_and_phi_deim_output_dofs_[cell];
    phi_output_dofs.clear();
    phinat_output_dofs.clear();
    mu_output_dofs.clear();
    phinat_mu_output_dofs.clear();
    for (const auto& dof : unique_output_dofs) {
      if (dof < size_phi_)
        phi_output_dofs.push_back(dof);
      else if (dof < 2 * size_phi_)
        phinat_output_dofs.push_back(dof);
      else
        mu_output_dofs.push_back(dof);
    }
    for (auto& dof : phinat_output_dofs)
      dof -= size_phi_;
    for (auto& dof : mu_output_dofs)
      dof -= 2 * size_phi_;
    for (const auto& dof : phinat_output_dofs)
      if (std::find(mu_output_dofs.begin(), mu_output_dofs.end(), dof) != mu_output_dofs.end())
        phinat_mu_output_dofs.push_back(dof);
    // get input dofs corresponding to output dofs
    auto& input_dofs = pfield_deim_input_dofs_[cell];
    input_dofs.clear();
    for (const auto& dof : unique_output_dofs) {
      const auto& new_input_dofs = pattern.inner(dof);
      input_dofs.insert(input_dofs.end(), new_input_dofs.begin(), new_input_dofs.end());
    }
    // sort and remove duplicate entries
    std::sort(input_dofs.begin(), input_dofs.end());
    input_dofs.erase(std::unique(input_dofs.begin(), input_dofs.end()), input_dofs.end());
    phinat_deim_input_dofs_begin_[cell] =
        std::lower_bound(input_dofs.begin(), input_dofs.end(), size_phi_) - input_dofs.begin();
    mu_deim_input_dofs_begin_[cell] =
        std::lower_bound(input_dofs.begin(), input_dofs.end(), 2 * size_phi_) - input_dofs.begin();
    // store all entities that contain an output dof
    const auto& mapper = phi_space_.mapper();
    DynamicVector<size_t> global_indices;
    pfield_deim_entities_[cell].clear();
    for (const auto& entity : Dune::elements(grid_view_)) {
      mapper.global_indices(entity, global_indices);
      maybe_add_entity(entity, global_indices, *pfield_deim_output_dofs_[cell], pfield_deim_entities_[cell], size_phi_);
    } // entities
  } // if (not already computed)
} // void compute_restricted_pfield_dofs(...)

//******************************************************************************************************************
//*********************************************** Apply operators **************************************************
//******************************************************************************************************************

// Applies stokes operator (applies the F if Stokes equation is F(y) = 0)
CellModelSolver::VectorType CellModelSolver::apply_stokes_operator(VectorType y, const bool /*restricted*/) const
{
  VectorType ret(size_u_ + size_p_, 0.);
  S_stokes_.mv(y, ret);
  ret -= stokes_rhs_vector_;
  return ret;
}

void CellModelSolver::assemble_nonlinear_part_of_ofield_residual(VectorType& residual,
                                                                 const size_t cell,
                                                                 const bool restricted)
{
  // nonlinear part
  VectorViewType res1_vec(residual, size_u_, 2 * size_u_);
  auto nonlinear_res_functional = make_vector_functional(u_space_, res1_vec);
  XT::Functions::GenericGridFunction<E, d, 1> nonlinear_res_pf(
      /*order = */ 3 * u_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_P(cell, element); },
      /*evaluate_func*/
      [cell, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate P, divP
        auto ret = this->eval_P(cell, x_local, param);
        ret *= factor * (ret * ret);
        return ret;
      });
  nonlinear_res_functional.append(LocalElementIntegralFunctional<E, d>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, d>(), nonlinear_res_pf)));
  if (!restricted)
    nonlinear_res_functional.assemble(use_tbb_);
  else
    nonlinear_res_functional.assemble_range(ofield_deim_entities_[cell]);
}

// Applies cell-th orientation field operator (applies F if the orientation field equation is F(y) = 0)
CellModelSolver::VectorType
CellModelSolver::apply_ofield_operator(const VectorType& y, const size_t cell, const bool restricted)
{
  const auto& output_dofs = *ofield_deim_output_dofs_[cell];
  const auto& unique_output_dofs = ofield_deim_unique_output_dofs_[cell];
  const auto& input_dofs = ofield_deim_input_dofs_[cell];
  auto& source = ofield_tmp_vec_;
  auto& residual = ofield_tmp_vec2_;
  // copy values to high-dimensional vector
  if (restricted)
    copy_ld_to_hd_vec(input_dofs, y, source);
  else
    source = y;
  // linear part
  ofield_jac_linear_op_.apply(source, residual);
  // subtract rhs
  const auto sub = sub_func<VectorType>(restricted);
  sub(residual, ofield_rhs_vector_, unique_output_dofs);
  // nonlinear part
  fill_tmp_ofield(cell, source, restricted);
  assemble_nonlinear_part_of_ofield_residual(residual, cell, restricted);
  if (restricted) {
    VectorType ret(output_dofs.size());
    for (size_t ii = 0; ii < output_dofs.size(); ++ii)
      ret[ii] = residual[output_dofs[ii]];
    return ret;
  } else {
    return residual;
  }
}

void CellModelSolver::update_ofield_parameters(const double Pa)
{
  // Pa may have been set to a new value already (via update_pfield_parameters)
  if (XT::Common::FloatCmp::ne(Pa, last_ofield_Pa_)) {
    std::cout << "Ofield params updated, old Pa = " << last_ofield_Pa_ << ", new Pa = " << Pa << std::endl;
    Pa_ = Pa;
    last_ofield_Pa_ = Pa_;
    C_ofield_elliptic_part_ *= 0.;
    MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, u_space_, u_space_, C_ofield_elliptic_part_);
    ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
    ofield_elliptic_op.assemble(use_tbb_);
  }
}

// Applies cell-th phase field operator (applies F if phase field equation is F(y) = 0)
CellModelSolver::VectorType
CellModelSolver::apply_pfield_operator(const VectorType& y, const size_t cell, const bool restricted)
{
  const auto& output_dofs = *pfield_deim_output_dofs_[cell];
  const auto& unique_output_dofs = pfield_deim_unique_output_dofs_[cell];
  const auto& input_dofs = pfield_deim_input_dofs_[cell];
  auto& source = pfield_tmp_vec_;
  auto& residual = pfield_tmp_vec2_;
  // copy values to high-dimensional vector
  if (restricted)
    copy_ld_to_hd_vec(input_dofs, y, source);
  else
    source = y;
  // linear part
  pfield_jac_linear_op_.apply(source, residual);
  // subtract rhs
  const auto sub = sub_func<VectorType>(restricted);
  sub(residual, pfield_rhs_vector_, unique_output_dofs);
  // nonlinear part
  fill_tmp_pfield(cell, source, restricted);
  assemble_nonlinear_part_of_pfield_residual(residual, cell, restricted);
  if (restricted) {
    VectorType ret(output_dofs.size());
    for (size_t ii = 0; ii < output_dofs.size(); ++ii)
      ret[ii] = residual[output_dofs[ii]];
    return ret;
  } else {
    return residual;
  }
}

void CellModelSolver::update_pfield_parameters(const double Be, const double Ca, const double Pa)
{
  if (XT::Common::FloatCmp::ne(Be, Be_) || XT::Common::FloatCmp::ne(Ca, Ca_)
      || XT::Common::FloatCmp::ne(Pa, last_pfield_Pa_)) {
    std::cout << "Pfield params updated, old (Be, Ca, Pa) = " << Be_ << ", " << Ca_ << ", " << last_pfield_Pa_
              << ", new = " << Be << ", " << Ca << ", " << Pa << std::endl;
    Be_ = Be;
    Ca_ = Ca;
    Pa_ = Pa;
    last_pfield_Pa_ = Pa_;
    XT::Common::Parameter param({{"gamma", {gamma_}}, {"epsilon", {epsilon_}}, {"Be", {Be_}}, {"Ca", {Ca_}}});
    pfield_jac_linear_op_.set_params(param);
    pfield_solver_.set_params(param);
  }
}

//******************************************************************************************************************
//******************************************* Apply inverse operators **********************************************
//******************************************************************************************************************

// Applies inverse stokes operator (solves F(y) = 0)
CellModelSolver::VectorType CellModelSolver::apply_inverse_stokes_operator() const
{
  // now solve the system
  // auto begin = std::chrono::steady_clock::now();
  EigenVectorType ret(size_u_ + size_p_);
  ret.backend() = stokes_solver_->solve(stokes_rhs_vector_.backend());
  // std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
  // std::cout << "Solving Stokes took: " << time.count() << " s!" << std::endl;

  // ensure int_\Omega p = 0 (TODO: remove (?), not necessary as p is not used anywhere)
  // auto p_integral = p_basis_integrated_vector_ * p_.dofs().vector();
  // auto p_correction = make_discrete_function<VectorType>(p_space_, "p_corr");
  // XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain_);
  // default_interpolation(const_p_integral_func, p_correction);
  // p_ -= p_correction;

  return XT::Common::convert_to<VectorType>(ret);
}

// Applies inverse orientation field operator (solves F(y) = 0)
// y_guess is the initial guess for the Newton iteration
CellModelSolver::VectorType CellModelSolver::apply_inverse_ofield_operator(const VectorType& y_guess, const size_t cell)
{
  if (linearize_) {
    return ofield_solver_.apply(ofield_rhs_vector_, cell);
  } else {

    // *********** Newton ******************************
    const auto tol = 1e-10;
    const auto max_iter = 200;
    const auto max_dampening_iter = 1000;

    auto l2_norm_P = l2_norm(grid_view_, P_[cell]);
    auto l2_norm_Pnat = l2_norm(grid_view_, Pnat_[cell]);

    // ********* compute residual *********
    auto begin = std::chrono::steady_clock::now();
    auto residual = apply_ofield_operator(y_guess, cell);
    auto res_norm = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    // std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

    size_t iter = 0;
    VectorType y_n = y_guess;
    VectorType y_n_plus_1 = y_guess;
    VectorType update;
    while (true) {
      if (res_norm < tol)
        break;

      // ********** assemble nonlinear part of S = Jacobian ***********
      begin = std::chrono::steady_clock::now();
      assemble_ofield_nonlinear_jacobian(y_n, cell);
      time = std::chrono::steady_clock::now() - begin;
      // std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

      // *********** solve system *************
      residual *= -1.;
      update = ofield_solver_.apply(residual, cell);

      DUNE_THROW_IF(iter >= max_iter,
                    Exceptions::operator_error,
                    "max iterations in ofield reached!\n|residual|_l2 = " << res_norm << ", param: "
                                                                          << "(" << Re_ << ", " << 1. / Fa_inv_ << ", "
                                                                          << xi_ << ")" << std::endl);

      // apply damping
      size_t k = 0;
      auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
      double lambda = 1;

      // revert jacobian back to linear part to correctly calculate linear part of residual
      // revert_ofield_jacobian_to_linear();

      // backtracking line search
      const double gamma = 0.001;
      while (candidate_res > (1 - gamma * lambda) * res_norm) {
        DUNE_THROW_IF(k >= max_dampening_iter,
                      Exceptions::operator_error,
                      "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                          << res_norm << "\nl = " << iter << "\n");
        y_n_plus_1 = y_n + update * lambda;
        residual = apply_ofield_operator(y_n_plus_1, cell);
        candidate_res = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
        // std::cout << "Candidate res: " << candidate_res << std::endl;
        lambda /= 2;
        k += 1;
      }
      y_n = y_n_plus_1;
      res_norm = candidate_res;
      // std::cout << "Current res: " << candidate_res << std::endl;
      iter += 1;
    } // while (true)
    return y_n;
  }
}

// Applies inverse phase field operator (solves F(y) = 0)
// y_guess is the initial guess for the Newton iteration
CellModelSolver::VectorType CellModelSolver::apply_inverse_pfield_operator(const VectorType& y_guess, const size_t cell)
{
  if (linearize_) {
    return pfield_solver_.apply(pfield_rhs_vector_, cell);
  } else {

    // *********** Newton ******************************
    const auto tol = 1e-10;
    const auto max_iter = 200;
    const auto max_dampening_iter = 1000;

    const auto l2_norm_phi = l2_norm(grid_view_, phi_[cell]);
    const auto l2_norm_phinat = l2_norm(grid_view_, phinat_[cell]);
    const auto l2_norm_mu = l2_norm(grid_view_, mu_[cell]);

    // ********* compute residual *********
    auto begin = std::chrono::steady_clock::now();
    auto residual = apply_pfield_operator(y_guess, cell, false);
    auto res_norm = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    // std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

    size_t iter = 0;
    VectorType x_n = y_guess;
    VectorType x_n_plus_1 = y_guess;
    VectorType update;
    while (true) {
      if (res_norm < tol)
        break;

      // ********** assemble nonlinear part of S = Jacobian ***********
      begin = std::chrono::steady_clock::now();
      assemble_pfield_nonlinear_jacobian(x_n, cell, false);
      time = std::chrono::steady_clock::now() - begin;
      // std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

      // *********** solve system *************
      residual *= -1.;
      update = pfield_solver_.apply(residual, cell);

      DUNE_THROW_IF(iter >= max_iter,
                    Exceptions::operator_error,
                    "max iterations in pfield reached!\n|residual|_l2 = " << res_norm << ", param: "
                                                                          << "(" << Re_ << ", " << 1. / Fa_inv_ << ", "
                                                                          << xi_ << ")" << std::endl);

      // apply damping
      size_t k = 0;
      auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
      double lambda = 1;

      // revert jacobian back to linear part to correctly calculate linear part of residual
      // revert_pfield_jacobian_to_linear();

      // backtracking line search
      const double gamma = 0.001;
      while (candidate_res > (1 - gamma * lambda) * res_norm) {
        DUNE_THROW_IF(k >= max_dampening_iter,
                      Exceptions::operator_error,
                      "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                          << res_norm << "\nl = " << iter << "\n");
        x_n_plus_1 = x_n + update * lambda;
        residual = apply_pfield_operator(x_n_plus_1, cell, false);
        candidate_res = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
        // std::cout << "Candidate res: " << candidate_res << std::endl;
        lambda /= 2;
        k += 1;
      }
      x_n = x_n_plus_1;
      res_norm = candidate_res;
      // std::cout << "Current res: " << candidate_res << std::endl;
      iter += 1;
    } // while (true)
    return x_n;
  }
}

//******************************************************************************************************************
//********************************************** Apply jacobians ***************************************************
//******************************************************************************************************************

void CellModelSolver::set_pfield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted)
{
  assemble_pfield_nonlinear_jacobian(source, cell, restricted);
}

// Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
// As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
CellModelSolver::VectorType
CellModelSolver::apply_pfield_jacobian(const VectorType& source, const size_t cell, const bool restricted)
{
  VectorType range(source.size(), 0.);
  const auto& output_dofs = *pfield_deim_output_dofs_[cell];
  const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
  VectorType& full_range = restricted ? pfield_tmp_vec_ : range;
  VectorViewType range_phi(full_range, 0, size_phi_);
  VectorViewType range_phinat(full_range, size_phi_, 2 * size_phi_);
  VectorViewType range_mu(full_range, 2 * size_phi_, 3 * size_phi_);
  const ConstVectorViewType source_phi(source, 0, size_phi_);
  const ConstVectorViewType source_phinat(source, size_phi_, 2 * size_phi_);
  const ConstVectorViewType source_mu(source, 2 * size_phi_, 3 * size_phi_);
  // linear part
  pfield_jac_linear_op_.apply(source, full_range);
  // nonlinear_part
  auto& tmp_vec = phi_tmp_vec_;
  const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
  const auto axpy = axpy_func<VectorViewType, VectorType>(restricted);
  const auto scal = scal_func<VectorType>(restricted);
  const auto add = add_func<VectorViewType, VectorType>(restricted);
  // apply missing parts of J (including the linear 1./Ca_ part)
  mv(M_nonlin_pfield_, source_mu, tmp_vec, phinat_output_dofs);
  axpy(range_phinat, 1. / (Be_ * std::pow(epsilon_, 2)), tmp_vec, phinat_output_dofs);
  mv(M_pfield_, source_mu, tmp_vec, phinat_output_dofs);
  axpy(range_phinat, 1. / Ca_, tmp_vec, phinat_output_dofs);
  // apply missing parts of A
  mv(M_nonlin_pfield_, source_phi, tmp_vec, mu_output_dofs);
  scal(tmp_vec, 1. / epsilon_, mu_output_dofs);
  add(range_mu, tmp_vec, mu_output_dofs);
  // apply G
  mv(G_pfield_, source_phi, tmp_vec, phinat_output_dofs);
  add(range_phinat, tmp_vec, phinat_output_dofs);

  if (restricted)
    for (size_t ii = 0; ii < output_dofs.size(); ++ii)
      range.set_entry(ii, full_range.get_entry(output_dofs[ii]));
  return range;
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_pfield_jacobian(const VectorType& rhs, const size_t cell)
{
  return pfield_solver_.apply(rhs, cell);
}

void CellModelSolver::set_ofield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted)
{
  assemble_ofield_nonlinear_jacobian(source, cell, restricted);
}

// Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
// As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
CellModelSolver::VectorType
CellModelSolver::apply_ofield_jacobian(const VectorType& source, const size_t cell, const bool restricted)
{
  VectorType range(source.size(), 0.);
  const auto& output_dofs = *ofield_deim_output_dofs_[cell];
  const auto& Pnat_output_dofs = Pnat_deim_output_dofs_[cell];
  VectorType& full_range = restricted ? ofield_tmp_vec_ : range;
  VectorViewType range_P(full_range, 0, size_u_);
  VectorViewType range_Pnat(full_range, size_u_, 2 * size_u_);
  const ConstVectorViewType source_P(source, 0, size_u_);
  const ConstVectorViewType source_Pnat(source, size_u_, 2 * size_u_);
  // linear part
  ofield_jac_linear_op_.apply(source, full_range);
  // nonlinear_part
  auto& tmp_vec = u_tmp_vec_;
  const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
  const auto add = add_func<VectorViewType, VectorType>(restricted);
  mv(C_ofield_nonlinear_part_, source_P, tmp_vec, Pnat_output_dofs);
  add(range_Pnat, tmp_vec, Pnat_output_dofs);
  if (restricted)
    for (size_t ii = 0; ii < output_dofs.size(); ++ii)
      range.set_entry(ii, full_range.get_entry(output_dofs[ii]));
  return range;
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_ofield_jacobian(const VectorType& rhs, const size_t cell)
{
  return ofield_solver_.apply(rhs, cell);
}

CellModelSolver::VectorType CellModelSolver::apply_stokes_jacobian(const VectorType& source, const bool /*restricted*/)
{
  VectorType range(source.size());
  S_stokes_.mv(source, range);
  return range;
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_stokes_jacobian(const VectorType& rhs)
{
  EigenVectorType rhs_eigen = XT::Common::convert_to<EigenVectorType>(rhs);
  EigenVectorType ret(stokes_solver_->solve(rhs_eigen.backend()));
  return XT::Common::convert_to<VectorType>(ret);
}

//******************************************************************************************************************
//**************************** Methods to assemble rhs, residuals and jacobians ************************************
//******************************************************************************************************************

// Computes stokes rhs using currently stored values of variables and stores in stokes_rhs_vector_
void CellModelSolver::assemble_stokes_rhs()
{
  auto f_functional = make_vector_functional(u_space_, stokes_f_vector_);
  // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
  f_functional.append(LocalElementIntegralFunctional<E, d>(
      /*order*/ [& u_space =
                     u_space_](const auto& test_basis,
                               const auto& param) { return 3 * u_space.max_polorder() + test_basis.order(param); },
      /*evaluate_func*/
      [this](const auto& test_basis,
             const DomainType& x_local,
             DynamicVector<R>& result,
             const XT::Common::Parameter& param) {
        const size_t sz = test_basis.size(param);
        if (result.size() < sz)
          result.resize(sz);
        std::fill(result.begin(), result.end(), 0.);
        using TestBasisType = typename LocalElementIntegralFunctional<E, d>::GenericIntegrand::LocalBasisType;
        thread_local std::vector<typename TestBasisType::RangeType> test_basis_values_;
        thread_local std::vector<typename TestBasisType::DerivativeRangeType> test_basis_grads_;
        test_basis.evaluate(x_local, test_basis_values_, param);
        test_basis.jacobians(x_local, test_basis_grads_, param);

        // evaluate P, Pnat, phi, phinat, \nabla P, \nabla Pnat, \nabla phi, div P, div Pnat and phi_tilde = (phi +
        // 1)/2 return type of the jacobians is a FieldMatrix<r, d>
        for (size_t kk = 0; kk < this->num_cells_; ++kk) {
          const auto P = this->eval_P(kk, x_local, param);
          const auto Pnat = this->eval_Pnat(kk, x_local, param);
          const auto phi = this->eval_phi(kk, x_local, param);
          const auto phinat = this->eval_phinat(kk, x_local, param);
          const auto grad_P = this->grad_P(kk, x_local, param);
          const auto grad_phi = this->grad_phi(kk, x_local, param);
          const auto phi_tilde = (phi + 1.) / 2.;

          // evaluate rhs terms
          const auto phinat_grad_phi = grad_phi * phinat;
          auto grad_P_T_times_Pnat = P;
          grad_P.mtv(Pnat, grad_P_T_times_Pnat);
          for (size_t ii = 0; ii < sz; ++ii) {
            for (size_t mm = 0; mm < d; ++mm) {
              result[ii] += (phinat_grad_phi[mm] + grad_P_T_times_Pnat[mm]) * test_basis_values_[ii][mm];
              for (size_t nn = 0; nn < d; ++nn)
                result[ii] += (-this->Fa_inv_ * phi_tilde * P[mm] * P[nn] - 0.5 * (this->xi_ + 1) * Pnat[mm] * P[nn]
                               - 0.5 * (this->xi_ - 1) * P[mm] * Pnat[nn])
                              * test_basis_grads_[ii][mm][nn];
            } // mm
          } // ii
        } // kk
      },
      /*post_bind_func*/
      [this](const E& element) {
        for (size_t kk = 0; kk < this->num_cells_; ++kk) {
          this->bind_phi(kk, element);
          this->bind_phinat(kk, element);
          this->bind_P(kk, element);
          this->bind_Pnat(kk, element);
        } // kk
      }));
  stokes_f_vector_ *= 0.;
  f_functional.assemble(use_tbb_);
  // apply dirichlet constraints for u.
  u_dirichlet_constraints_.apply(stokes_f_vector_);
}

// Computes orientation field rhs using currently stored values of variables and stores in ofield_rhs_vector_
void CellModelSolver::assemble_ofield_rhs(const double /*dt*/, const size_t cell)
{
  auto g_functional = make_vector_functional(u_space_, ofield_g_vector_);
  ofield_g_vector_ *= 0.;
  XT::Functions::GenericGridFunction<E, d> g(
      /*order = */ 3 * u_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) {
        this->bind_phi(cell, element);
        if (this->linearize_) {
          this->bind_P(cell, element);
        }
      },
      /*evaluate_func*/
      [cell, factor1 = beta_ / Pa_, factor2 = -2. * c_1_ / Pa_, this](const DomainType& x_local,
                                                                      const XT::Common::Parameter& param) {
        // evaluate rhs terms
        const auto grad_phi = this->grad_phi(cell, x_local, param);
        auto ret = grad_phi;
        ret *= factor1;
        if (linearize_) {
          const auto P_n = this->eval_P(cell, x_local, param);
          auto ret2 = P_n;
          ret2 *= factor2 * (P_n * P_n);
          ret += ret2;
        }
        return ret;
      });
  g_functional.append(LocalElementIntegralFunctional<E, d>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, d>(), g)));
  g_functional.assemble(use_tbb_);
  M_ofield_.mv(P_[cell].dofs().vector(), ofield_f_vector_);
}

// Computes phase field rhs using currently stored values of variables and stores in pfield_rhs_vector_
void CellModelSolver::assemble_pfield_rhs(const double /*dt*/, const size_t cell, const bool restricted)
{
  const auto& phi_output_dofs = phi_deim_output_dofs_[cell];
  const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  auto f_functional = make_vector_functional(phi_space_, pfield_f_vector_);
  auto h_functional = make_vector_functional(phi_space_, pfield_h_vector_);
  // calculate f
  if (linearize_)
    pfield_f_vector_ *= 0.;
  XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
      /*order = */ 3 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_phi(cell, element); },
      /*evaluate_func*/
      [cell, two_epsilon_inv = 2. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate phi_
        const R phi_n = this->eval_phi(cell, x_local, param);
        return two_epsilon_inv * std::pow(phi_n, 3);
      });
  f_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), f_pf)));
  if (linearize_)
    h_functional.append(f_functional);
  // handling of dirichlet boundary conditions is not necessary as these are applied during computaton of the residual

  // calculate g
  const auto mv = mv_func<VectorViewType>(restricted);
  mv(M_pfield_, phi_[cell].dofs().vector(), pfield_g_vector_, phi_output_dofs);
  for (const auto& DoF : phi_dirichlet_constraints_.dirichlet_DoFs())
    pfield_g_vector_[DoF] = -1.;

  // calculate h
  const auto scal = scal_func<VectorViewType>(restricted);
  scal(pfield_h_vector_, 0., phinat_output_dofs);
  XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
      /*order = */ 3 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) {
        this->bind_P(cell, element);
        if (this->linearize_) {
          this->bind_phi(cell, element);
          this->bind_mu(cell, element);
        }
      },
      /*evaluate_func*/
      [cell, factor0 = 6. / (Be_ * std::pow(epsilon_, 2)), factor1 = -c_1_ / (2. * Pa_), factor2 = -beta_ / Pa_, this](
          const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate P, divP
        const auto Pn = this->eval_P(cell, x_local, param);
        const auto grad_P = this->grad_P(cell, x_local, param);
        R div_P(0.);
        for (size_t ii = 0; ii < d; ++ii)
          div_P += grad_P[ii][ii];
        auto ret = factor1 * (Pn * Pn) + factor2 * div_P;
        if (this->linearize_) {
          const auto phi_n = this->eval_phi(cell, x_local, param);
          const auto mu_n = this->eval_mu(cell, x_local, param);
          ret += factor0 * std::pow(phi_n, 2) * mu_n;
        }
        return ret;
      });
  h_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(1.), h_pf)));
  // assemble rhs
  if (!restricted)
    h_functional.assemble(use_tbb_);
  else
    h_functional.assemble_range(pfield_deim_entities_[cell]);
}

// assembles linear part of orientation field jacobian and stores in S_ofield_
void CellModelSolver::assemble_ofield_linear_jacobian(const double dt, const size_t cell)
{
  // calculate A
  // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
  XT::Functions::GenericGridFunction<E, d, d> Omega_minus_xi_D_transposed(
      /*order = */ std::max(u_space_.max_polorder() - 1, 0),
      /*post_bind_func*/
      [this](const E& element) { this->bind_u(element); },
      /*evaluate_func*/
      [this](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate \nabla u
        auto grad_u = this->grad_u(x_local, param);
        auto grad_u_T = grad_u;
        grad_u_T.transpose();
        auto& ret = grad_u;
        ret *= (1. - this->xi_) / 2.;
        grad_u_T *= (1. + this->xi_) / 2.;
        ret -= grad_u_T;
        return ret;
      });
  A_ofield_.set_to_zero();
  A_ofield_op_->clear();
  A_ofield_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(Omega_minus_xi_D_transposed)));
  A_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(u_)));
  A_ofield_op_->assemble(use_tbb_);

  // calculate linear part S_10 = C
  C_ofield_linear_part_op_->clear();
  C_ofield_linear_part_ = C_ofield_elliptic_part_;
  XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_inv_phi(
      /*order = */ u_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_phi(cell, element); },
      /*evaluate_func*/
      [cell, factor = c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi = this->eval_phi(cell, x_local, param);
        return factor * phi;
      });
  C_ofield_linear_part_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(c1_Pa_inv_phi)));
  C_ofield_linear_part_op_->assemble(use_tbb_);
  if (ofield_solver_.is_schur_solver()) {
    S_schur_ofield_linear_part_.backend() = M_ofield_.backend();
    S_schur_ofield_linear_part_.axpy(dt, A_ofield_);
    S_schur_ofield_linear_part_.axpy(-dt / kappa_, C_ofield_linear_part_);
  }

  // nonlinear part is equal to linearized part in first iteration
  if (linearize_)
    assemble_ofield_nonlinear_jacobian(ofield_vec(cell), cell);
}

// assembles nonlinear part of orientation field jacobian and adds to S_ofield_
// if assemble_ofield_linear_jacobian has been called first, S_ofield now contains the whole orientation field
// jacobian
void CellModelSolver::assemble_ofield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  fill_tmp_ofield(cell, y, restricted);
  assemble_C_ofield_nonlinear_part(cell, restricted);
  if (ofield_solver_.is_schur_solver()) {
    ofield_solver_.schur_matrix() = S_schur_ofield_linear_part_;
    ofield_solver_.schur_matrix().axpy(-dt_ / kappa_, C_ofield_nonlinear_part_);
  }
  ofield_solver_.prepare(dt_, cell, restricted);
}

void CellModelSolver::assemble_C_ofield_nonlinear_part(const size_t cell, const bool restricted)
{
  XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_P2(
      /*order = */ 2. * u_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_P(cell, element); },
      /*evaluate_func*/
      [cell, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto P_n = this->eval_P(cell, x_local, param);
        return factor * P_n.two_norm2();
      });
  set_mat_to_zero(C_ofield_nonlinear_part_, restricted, ofield_submatrix_pattern_, Pnat_deim_output_dofs_[cell]);
  C_ofield_nonlinear_part_op_->clear();
  C_ofield_nonlinear_part_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(c1_Pa_P2)));
  C_ofield_nonlinear_part_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementOtimesMatrixIntegrand<E, d>(P_tmp_[cell], -2. * c_1_ / Pa_)));
  if (!restricted)
    C_ofield_nonlinear_part_op_->assemble(use_tbb_);
  else
    C_ofield_nonlinear_part_op_->assemble_range(ofield_deim_entities_[cell]);
}

// assembles linear part of phase field jacobian
void CellModelSolver::assemble_pfield_linear_jacobian(const double /*dt*/, const size_t cell, const bool restricted)
{
  // assemble matrix S_{00} = M + dt D
  assemble_D_pfield(cell, restricted);
  // nonlinear part is equal to linearized part in first iteration
  if (linearize_)
    assemble_pfield_nonlinear_jacobian(pfield_vec(cell), cell, restricted);
}

void CellModelSolver::set_mat_to_zero(MatrixType& mat,
                                      const bool restricted,
                                      const XT::LA::SparsityPatternDefault& pattern,
                                      const std::vector<size_t>& rows)
{
  if (!restricted) {
    mat.set_to_zero();
  } else {
    for (const auto& row : rows)
      for (const auto& col : pattern.inner(row))
        mat.set_entry(row, col, 0.);
  }
}

void CellModelSolver::assemble_D_pfield(const size_t cell, const bool restricted)
{
  const auto& phi_output_dofs = phi_deim_output_dofs_[cell];
  set_mat_to_zero(D_pfield_, restricted, pfield_submatrix_pattern_, phi_output_dofs);
  XT::Functions::GenericGridFunction<E, d, 1> minus_u(
      /*order = */ u_space_.max_polorder(),
      /*post_bind_func*/
      [this](const E& element) { this->bind_u(element); },
      /*evaluate_func*/
      [this](const DomainType& x_local, const XT::Common::Parameter& param) {
        auto ret = this->eval_u(x_local, param);
        ret *= -1;
        return ret;
      });
  D_pfield_op_->clear();
  D_pfield_op_->append(
      LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u)));
  if (!restricted)
    D_pfield_op_->assemble(use_tbb_);
  else
    D_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
}

// assembles nonlinear part of phase field jacobian
void CellModelSolver::assemble_pfield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  fill_tmp_pfield(cell, y, restricted);
  assemble_M_nonlin_pfield(cell, restricted);
  assemble_G_pfield(cell, restricted);
  pfield_solver_.prepare(dt_, cell, restricted);
}

// stores matrix with entries \int (3 phi^2 - 1) varphi_i varphi_j in M_nonlin_pfield_
void CellModelSolver::assemble_M_nonlin_pfield(const size_t cell, const bool restricted)
{
  const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  const auto& mu_output_dofs = mu_deim_output_dofs_[cell];
  set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, mu_output_dofs);
  if (restricted)
    set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, phinat_output_dofs);
  XT::Functions::GenericGridFunction<E, 1, 1> M_nonlin_prefactor(
      /*order = */ 2 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_phi(cell, element); },
      /*evaluate_func*/
      [cell, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        const R phi_n = this->eval_phi(cell, x_local, param);
        return (3. * phi_n * phi_n - 1.);
      });
  M_nonlin_pfield_op_->clear();
  M_nonlin_pfield_op_->append(
      LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(M_nonlin_prefactor)));
  if (!restricted)
    M_nonlin_pfield_op_->assemble(use_tbb_);
  else
    M_nonlin_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
}

// stores nonlinear part of block G of the phase field jacobian matrix in G_pfield_nonlinear_part_
void CellModelSolver::assemble_G_pfield(const size_t cell, const bool restricted)
{
  const auto& phinat_output_dofs = phinat_deim_output_dofs_[cell];
  set_mat_to_zero(G_pfield_, restricted, pfield_submatrix_pattern_, phinat_output_dofs);
  XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
      /*order = */ 2 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) {
        this->bind_phi(cell, element);
        this->bind_mu(cell, element);
        if (this->num_cells_ > 1) {
          for (size_t kk = 0; kk < this->num_cells_; ++kk)
            this->bind_phi(kk, element);
        }
      },
      /*evaluate_func*/
      [cell, In_inv = 1. / In_, eps_inv = 1. / epsilon_, six_inv_Be_eps2 = 6. / (Be_ * std::pow(epsilon_, 2)), this](
          const DomainType& x_local, const XT::Common::Parameter& param) {
        const R phi_n = this->eval_phi(cell, x_local, param);
        const R mu_n = this->eval_mu(cell, x_local, param);
        auto ret = six_inv_Be_eps2 * phi_n * mu_n;
        if (this->num_cells_ > 1) {
          R wsum = 0.;
          R Bsum = 0.;
          for (size_t kk = 0; kk < this->num_cells_; ++kk) {
            if (kk != cell) {
              wsum += this->w_func(kk, x_local, param);
              Bsum += this->B_func(kk, x_local, param);
            }
          } // kk
          ret += In_inv * 4 * eps_inv * (3. * std::pow(phi_n, 2) - 1) * wsum;
          auto w_twoprime = 0;
          if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.)) {
            const auto ln = std::log((1 + phi_n) / (1 - phi_n));
            const auto ln2 = std::pow(ln, 2);
            w_twoprime = 4 * std::exp(-0.5 * ln2) * (ln2 - phi_n * ln - 1) / (std::pow(std::pow(phi_n, 2) - 1, 2));
          }
          ret += In_inv * w_twoprime * Bsum;
        } // num_cells > 1
        return ret;
      });
  G_pfield_op_->clear();
  G_pfield_op_->append(
      LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(G_prefactor)));
  if (!restricted)
    G_pfield_op_->assemble(use_tbb_);
  else
    G_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
}

// assembles nonlinear part of phasefield residual and adds to residual
void CellModelSolver::assemble_nonlinear_part_of_pfield_residual(VectorType& residual,
                                                                 const size_t cell,
                                                                 const bool restricted)
{
  VectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
  VectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);
  const auto res1 = make_discrete_function(phi_space_, res1_vec);
  const auto res2 = make_discrete_function(phi_space_, res2_vec);
  auto nonlinear_res1_functional = make_vector_functional(phi_space_, res1_vec);
  auto nonlinear_res2_functional = make_vector_functional(phi_space_, res2_vec);
  XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf1(
      /*order = */ 3 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) {
        this->bind_phi(cell, element);
        this->bind_mu(cell, element);
        if (this->num_cells_ > 1) {
          for (size_t kk = 0; kk < this->num_cells_; ++kk)
            this->bind_phi(kk, element);
        }
      },
      /*evaluate_func*/
      [cell,
       Ca_inv = 1. / Ca_,
       In_inv = 1. / In_,
       eps_inv = 1. / epsilon_,
       num_cells = num_cells_,
       inv_Be_eps2 = 1. / (Be_ * std::pow(epsilon_, 2)),
       this](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate P, divP
        const auto phi_n = this->eval_phi(cell, x_local, param);
        const auto mu_n = this->eval_mu(cell, x_local, param);
        auto ret = (Ca_inv + inv_Be_eps2 * (3. * phi_n * phi_n - 1)) * mu_n;
        if (num_cells > 1) {
          R wsum = 0.;
          R Bsum = 0.;
          for (size_t kk = 0; kk < num_cells; ++kk) {
            if (kk != cell) {
              wsum += this->w_func(kk, x_local, param);
              Bsum += this->B_func(kk, x_local, param);
            }
          } // kk
          ret += In_inv * 4 * eps_inv * (std::pow(phi_n, 3) - phi_n) * wsum;
          auto w_prime = 0;
          if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.)) {
            const auto ln = std::log((1 + phi_n) / (1 - phi_n));
            const auto ln2 = std::pow(ln, 2);
            w_prime = 2 * std::exp(-0.5 * ln2) * ln / (std::pow(phi_n, 2) - 1);
          }
          ret += In_inv * w_prime * Bsum;
        } // num_cells > 1
        return ret;
      });
  XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf2(
      /*order = */ 3 * phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_phi(cell, element); },
      /*evaluate_func*/
      [cell, inv_eps = 1. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        // evaluate P, divP
        const auto phi_n = this->eval_phi(cell, x_local, param);
        return inv_eps * (phi_n * phi_n - 1) * phi_n;
      });
  nonlinear_res1_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), nonlinear_res_pf1)));
  nonlinear_res2_functional.append(LocalElementIntegralFunctional<E, 1>(
      local_binary_to_unary_element_integrand(LocalElementProductScalarWeightIntegrand<E, 1>(), nonlinear_res_pf2)));
  nonlinear_res1_functional.append(nonlinear_res2_functional);
  if (!restricted)
    nonlinear_res1_functional.assemble(use_tbb_);
  else
    nonlinear_res1_functional.walk_range(pfield_deim_entities_[cell]);
  // high-dimensional operation, TODO: replace
  // phi_dirichlet_constraints_.apply(res2_vec);
}

//******************************************************************************************************************
//******************************************* DEIM related methods *************************************************
//******************************************************************************************************************

// Dofs needed for evaluation of output_dofs provided in
std::vector<size_t> CellModelSolver::pfield_deim_input_dofs(const size_t cell) const
{
  return pfield_deim_input_dofs_[cell];
}

size_t CellModelSolver::pfield_deim_input_dofs_size(const size_t cell) const
{
  return pfield_deim_input_dofs_[cell].size();
}

// Dofs needed for evaluation of output_dofs provided in
std::vector<size_t> CellModelSolver::ofield_deim_input_dofs(const size_t cell) const
{
  return ofield_deim_input_dofs_[cell];
}

// private:
//******************************************************************************************************************
//************ The following methods all bind or evaluate the respective temporary discrete function ***************
//******************************************************************************************************************

void CellModelSolver::bind_u(const E& element) const
{
  auto& u_local = **u_tmp_local_;
  if (!u_local)
    u_local = u_tmp_.local_function();
  u_local->bind(element);
}

CellModelSolver::DomainRetType CellModelSolver::eval_u(const DomainType& x_local,
                                                       const XT::Common::Parameter& param) const
{
  return (**u_tmp_local_)->evaluate(x_local, param);
}

CellModelSolver::JacobianRetType CellModelSolver::grad_u(const DomainType& x_local,
                                                         const XT::Common::Parameter& param) const
{
  return (**u_tmp_local_)->jacobian(x_local, param);
}

void CellModelSolver::bind_P(const size_t cell, const E& element) const
{
  auto& P_local_cell = (**P_tmp_local_)[cell];
  if (!P_local_cell)
    P_local_cell = P_tmp_[cell].local_function();
  P_local_cell->bind(element);
}

CellModelSolver::DomainRetType
CellModelSolver::eval_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
{
  auto& P_local_cell = (**P_tmp_local_)[cell];
  return P_local_cell->evaluate(x_local, param);
}

CellModelSolver::JacobianRetType
CellModelSolver::grad_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
{
  auto& P_local_cell = (**P_tmp_local_)[cell];
  return P_local_cell->jacobian(x_local, param);
}

void CellModelSolver::bind_Pnat(const size_t cell, const E& element) const
{
  auto& Pnat_local_cell = (**Pnat_tmp_local_)[cell];
  if (!Pnat_local_cell)
    Pnat_local_cell = Pnat_tmp_[cell].local_function();
  Pnat_local_cell->bind(element);
}

CellModelSolver::DomainRetType
CellModelSolver::eval_Pnat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const
{
  auto& Pnat_local_cell = (**Pnat_tmp_local_)[cell];
  return Pnat_local_cell->evaluate(x_local, param);
}

void CellModelSolver::bind_phi(const size_t cell, const E& element) const
{
  auto& phi_local_cell = (**phi_tmp_local_)[cell];
  if (!phi_local_cell)
    phi_local_cell = phi_tmp_[cell].local_function();
  phi_local_cell->bind(element);
}

CellModelSolver::R
CellModelSolver::eval_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
{
  auto& phi_local_cell = (**phi_tmp_local_)[cell];
  return phi_local_cell->evaluate(x_local, param)[0];
}

CellModelSolver::DomainRetType
CellModelSolver::grad_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
{
  auto& phi_local_cell = (**phi_tmp_local_)[cell];
  return phi_local_cell->jacobian(x_local, param)[0];
}

void CellModelSolver::bind_phinat(const size_t cell, const E& element) const
{
  auto& phinat_local_cell = (**phinat_tmp_local_)[cell];
  if (!phinat_local_cell)
    phinat_local_cell = phinat_tmp_[cell].local_function();
  phinat_local_cell->bind(element);
}

CellModelSolver::R
CellModelSolver::eval_phinat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
{
  auto& phinat_local_cell = (**phinat_tmp_local_)[cell];
  return phinat_local_cell->evaluate(x_local, param)[0];
}

void CellModelSolver::bind_mu(const size_t cell, const E& element) const
{
  auto& mu_local_cell = (**mu_tmp_local_)[cell];
  if (!mu_local_cell)
    mu_local_cell = mu_tmp_[cell].local_function();
  mu_local_cell->bind(element);
}

CellModelSolver::R
CellModelSolver::eval_mu(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param)
{
  auto& mu_local_cell = (**mu_tmp_local_)[cell];
  return mu_local_cell->evaluate(x_local, param)[0];
}

//******************************************************************************************************************
//****************  Linear algebra operations acting only on parts of the given matrices and vectors ***************
//******************************************************************************************************************

// Copies low-dimensional vec to given entries of high-dimensional vec.
void CellModelSolver::copy_ld_to_hd_vec(const std::vector<size_t> dofs, const VectorType& ld_vec, VectorType& hd_vec)
{
  assert(ld_vec.size() == dofs.size());
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    hd_vec[dofs[ii]] = ld_vec[ii];
}

//***************************************************************************************************************
//*********************************************  Helper methods  ************************************************
//***************************************************************************************************************

// get lower left of computational domain from testcase name
XT::Common::FieldVector<CellModelSolver::R, CellModelSolver::d>
CellModelSolver::get_lower_left(const std::string& testcase)
{
  if (testcase == "single_cell")
    return {{0., 0.}};
  else if (testcase == "two_cells")
    return {{0., 0.}};
  else
    DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return FieldVector<R, d>();
}

// get upper right of computational domain from testcase name
XT::Common::FieldVector<CellModelSolver::R, CellModelSolver::d>
CellModelSolver::get_upper_right(const std::string& testcase)
{
  if (testcase == "single_cell")
    return {{160., 40.}};
  else if (testcase == "two_cells")
    return {{50., 50.}};
  else
    DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return FieldVector<R, d>();
}

// get directions in which domain is periodic from testcase name
std::string CellModelSolver::get_periodic_directions(const std::string& testcase)
{
  if (testcase == "single_cell")
    return "01";
  else if (testcase == "two_cells")
    return "00";
  else
    DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return "";
}

// get number of cells from testcase name
size_t CellModelSolver::get_num_cells(const std::string& testcase)
{
  if (testcase == "single_cell")
    return 1;
  else if (testcase == "two_cells")
    return 2;
  else
    DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return 0;
}

// creates sparsity pattern of stokes system matrix
XT::LA::SparsityPatternDefault CellModelSolver::create_stokes_pattern(const SpaceInterface<PGV, d, 1, R>& u_space,
                                                                      const SpaceInterface<PGV, 1, 1, R>& p_space)
{
  const auto pattern_A = make_element_sparsity_pattern(u_space, u_space, u_space.grid_view());
  const auto pattern_B = make_element_sparsity_pattern(u_space, p_space, u_space.grid_view());
  const auto m = u_space.mapper().size();
  const auto n = p_space.mapper().size();
  XT::LA::SparsityPatternDefault pattern(m + n);
  for (size_t ii = 0; ii < m; ++ii) {
    for (const auto& jj : pattern_A.inner(ii))
      pattern.insert(ii, jj);
    for (const auto& jj : pattern_B.inner(ii)) {
      pattern.insert(ii, m + jj);
      pattern.insert(m + jj, ii);
    }
  }
  pattern.insert(m, m);
  pattern.sort();
  return pattern;
}

// appends entity to input_entities if one of its global_indices is in output_dofs
void CellModelSolver::maybe_add_entity(const E& entity,
                                       const DynamicVector<size_t>& global_indices,
                                       const std::vector<size_t>& output_dofs,
                                       std::vector<E>& input_entities,
                                       const size_t subvector_size) const
{
  for (const auto& output_dof : output_dofs) {
    const size_t dof = output_dof % subvector_size;
    for (size_t jj = 0; jj < global_indices.size(); ++jj) {
      if (global_indices[jj] == dof) {
        input_entities.push_back(entity);
        return;
      }
    } // jj
  } // dof
}

// sets temporary orientation field discrete functions to source values
void CellModelSolver::fill_tmp_ofield(const size_t cell, const VectorType& source, const bool restricted) const
{
  if (!restricted) {
    ConstVectorViewType P_vec(source, 0, size_u_);
    P_tmp_[cell].dofs().vector() = P_vec;
  } else {
    const auto& input_dofs = ofield_deim_input_dofs_[cell];
    for (size_t ii = 0; ii < Pnat_deim_input_dofs_begin_[cell]; ++ii)
      P_tmp_[cell].dofs().vector().set_entry(input_dofs[ii], source[input_dofs[ii]]);
  }
}

// sets temporary phase field discrete functions to source values
void CellModelSolver::fill_tmp_pfield(const size_t cell, const VectorType& source, const bool restricted) const
{
  if (!restricted) {
    ConstVectorViewType phi_vec(source, 0, size_phi_);
    ConstVectorViewType mu_vec(source, 2 * size_phi_, 3 * size_phi_);
    phi_tmp_[cell].dofs().vector() = phi_vec;
    mu_tmp_[cell].dofs().vector() = mu_vec;
  } else {
    const auto& input_dofs = pfield_deim_input_dofs_[cell];
    for (size_t ii = 0; ii < phinat_deim_input_dofs_begin_[cell]; ++ii)
      phi_tmp_[cell].dofs().vector().set_entry(input_dofs[ii], source[input_dofs[ii]]);
    for (size_t ii = mu_deim_input_dofs_begin_[cell]; ii < input_dofs.size(); ++ii)
      mu_tmp_[cell].dofs().vector().set_entry(input_dofs[ii] - 2 * size_phi_, source[input_dofs[ii]]);
  }
}

// error norm used in orientation field Newton iteration
// TODO: use appropriate norm
double CellModelSolver::ofield_residual_norm(const VectorType& residual, double l2_ref_P, double l2_ref_Pnat) const
{
  l2_ref_P = l2_ref_P < 1. ? 1. : l2_ref_P;
  l2_ref_Pnat = l2_ref_Pnat < 1. ? 1. : l2_ref_Pnat;
  ConstVectorViewType res0_vec(residual, 0, size_u_);
  ConstVectorViewType res1_vec(residual, size_u_, 2 * size_u_);
  const auto res0 = make_discrete_function(u_space_, res0_vec);
  const auto res1 = make_discrete_function(u_space_, res1_vec);
  return l2_norm(grid_view_, res0) / l2_ref_P + l2_norm(grid_view_, res1) / l2_ref_Pnat;
}

// error norm used in phase field Newton iteration
// TODO: use appropriate norm
double CellModelSolver::pfield_residual_norm(const VectorType& residual,
                                             double l2_ref_phi,
                                             double l2_ref_phinat,
                                             double l2_ref_mu) const
{
  l2_ref_phi = l2_ref_phi < 1. ? 1. : l2_ref_phi;
  l2_ref_phinat = l2_ref_phinat < 1. ? 1. : l2_ref_phinat;
  l2_ref_mu = l2_ref_mu < 1. ? 1. : l2_ref_mu;
  ConstVectorViewType res0_vec(residual, 0, size_phi_);
  ConstVectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
  ConstVectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);
  const auto res0 = make_discrete_function(phi_space_, res0_vec);
  const auto res1 = make_discrete_function(phi_space_, res1_vec);
  const auto res2 = make_discrete_function(phi_space_, res2_vec);
  return l2_norm(grid_view_, res0) / l2_ref_phi + l2_norm(grid_view_, res1) / l2_ref_phinat
         + l2_norm(grid_view_, res2) / l2_ref_mu;
}

CellModelSolver::R
CellModelSolver::B_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param)
{
  const R phi_n = eval_phi(kk, x_local, param);
  return 1 / epsilon_ * std::pow(std::pow(phi_n, 2) - 1, 2);
}

CellModelSolver::R
CellModelSolver::w_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param)
{
  const R phi_n = eval_phi(kk, x_local, param);
  if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
    return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
  else
    return 0.;
}