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

#include <dune/xt/common/numeric.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/filters.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/information.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>
#include <dune/xt/functions/visualization.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/symmetrized-laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/gradient-value.hh>
#include <dune/gdt/operators/bilinear-form.hh>
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
                                 const double dt,
                                 const unsigned int num_elements_x,
                                 const unsigned int num_elements_y,
                                 const int pol_order,
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
                                 const int inner_verbose)
  : lower_left_(get_lower_left(testcase))
  , upper_right_(get_upper_right(testcase))
  , t_end_(t_end)
  , dt_(dt)
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
  , In_(In)
  , vol_domain_((upper_right_[0] - lower_left_[0]) * (upper_right_[1] - lower_left_[1]))
  , num_cells_(get_num_cells(testcase))
  // do a global refine once, this makes simplicial grids look more symmetric
  , grid_(XT::Grid::make_cube_grid<G>(lower_left_, upper_right_, {num_elements_x, num_elements_y}, 1))
  , nonperiodic_grid_view_(grid_.leaf_view())
  , grid_view_(nonperiodic_grid_view_, std::bitset<d>(get_periodic_directions(testcase)))
  , dx_(XT::Grid::Dimensions<PGV>(grid_view_).entity_width.max())
  , epsilon_(epsilon < DXTC_CONFIG_GET("ratio_epsilon_dx", 3) * dx_ ? DXTC_CONFIG_GET("ratio_epsilon_dx", 3) * dx_
                                                                    : epsilon)
  , u_space_(make_continuous_lagrange_space<d>(grid_view_, 2))
  , P_space_(make_continuous_lagrange_space<d>(grid_view_, pol_order))
  , p_space_(make_continuous_lagrange_space<1>(grid_view_, 1))
  , phi_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order))
  , size_u_(u_space_.mapper().size())
  , size_P_(P_space_.mapper().size())
  , size_p_(p_space_.mapper().size())
  , size_phi_(phi_space_.mapper().size())
  , num_mutexes_u_(use_tbb ? 0 : size_u_ / 100)
  , num_mutexes_ofield_(use_tbb ? 0 : size_P_ / 100)
  , num_mutexes_pfield_(use_tbb ? 0 : size_phi_ / 100)
  , stokes_vector_(size_u_ + size_p_, 0., 0)
  , ofield_vectors_(num_cells_, VectorType(2 * size_P_, 0., 0))
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
  , stokes_tmp_vec_(size_u_ + size_p_)
  , stokes_tmp_vec2_(size_u_ + size_p_)
  , u_dirichlet_constraints_(make_dirichlet_constraints(u_space_, boundary_info_))
  , ofield_submatrix_pattern_(make_element_sparsity_pattern(P_space_, P_space_, grid_view_))
  , M_ofield_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , A_ofield_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_elliptic_part_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_linear_part_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_nonlinear_part_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , S_schur_ofield_linear_part_(size_P_, size_P_, ofield_submatrix_pattern_, 0)
  , M_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, P_space_, P_space_, M_ofield_))
  , A_ofield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, P_space_, P_space_, A_ofield_))
  , C_ofield_linear_part_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, P_space_, P_space_, C_ofield_linear_part_))
  , C_ofield_nonlinear_part_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, P_space_, P_space_, C_ofield_nonlinear_part_))
  , ofield_jac_linear_op_(M_ofield_, A_ofield_, C_ofield_linear_part_, dt_, kappa_, this)
  , ofield_solver_(dt_,
                   kappa_,
                   M_ofield_,
                   A_ofield_,
                   C_ofield_linear_part_,
                   C_ofield_nonlinear_part_,
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
  , ofield_rhs_vector_(2 * size_P_, 0., num_mutexes_ofield_)
  , ofield_f_vector_(ofield_rhs_vector_, 0, size_P_)
  , ofield_g_vector_(ofield_rhs_vector_, size_P_, 2 * size_P_)
  , stokes_solver_(std::make_shared<LUSolverType>())
  , ofield_tmp_vec_(2 * size_P_, 0., 0)
  , ofield_tmp_vec2_(2 * size_P_, 0., 0)
  , ofield_deim_source_dofs_(num_cells_)
  , Pnat_deim_source_dofs_begin_(num_cells_)
  , ofield_deim_range_dofs_(num_cells_)
  , ofield_deim_unique_range_dofs_(num_cells_)
  , P_deim_range_dofs_(num_cells_)
  , Pnat_deim_range_dofs_(num_cells_)
  , ofield_deim_entities_(num_cells_)
  , pfield_submatrix_pattern_(make_element_sparsity_pattern(phi_space_, phi_space_, grid_view_))
  , M_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , B_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , M_ell_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , M_nonlin_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , Dphi_f_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , B_pfield_op_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, B_pfield_))
  , Dphi_f_pfield_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, Dphi_f_pfield_))
  , M_nonlin_pfield_op_(
        std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(grid_view_, phi_space_, phi_space_, M_nonlin_pfield_))
  , pfield_jac_linear_op_(M_pfield_, B_pfield_, M_ell_pfield_, dt_, gamma_, epsilon_, Be_, this)
  , pfield_solver_(dt_,
                   gamma_,
                   epsilon_,
                   Be_,
                   Ca_,
                   M_pfield_,
                   M_ell_pfield_,
                   B_pfield_,
                   Dphi_f_pfield_,
                   M_nonlin_pfield_,
                   pfield_solver_type,
                   pfield_mass_matrix_solver_type,
                   pfield_submatrix_pattern_,
                   num_cells_,
                   outer_reduction,
                   outer_restart,
                   outer_verbose,
                   inner_reduction,
                   inner_maxit,
                   inner_verbose)
  , pfield_rhs_vector_(3 * size_phi_, 0., num_mutexes_pfield_)
  , pfield_r0_vector_(pfield_rhs_vector_, 0, size_phi_)
  , pfield_r1_vector_(pfield_rhs_vector_, size_phi_, 2 * size_phi_)
  , pfield_r2_vector_(pfield_rhs_vector_, 2 * size_phi_, 3 * size_phi_)
  , pfield_deim_source_dofs_(num_cells_)
  , phinat_deim_source_dofs_begin_(num_cells_)
  , mu_deim_source_dofs_begin_(num_cells_)
  , pfield_deim_range_dofs_(num_cells_)
  , pfield_deim_unique_range_dofs_(num_cells_)
  , phi_deim_range_dofs_(num_cells_)
  , phinat_deim_range_dofs_(num_cells_)
  , mu_deim_range_dofs_(num_cells_)
  , both_mu_and_phi_deim_range_dofs_(num_cells_)
  , pfield_deim_entities_(num_cells_)
  , pfield_tmp_vec_(3 * size_phi_, 0., 0)
  , pfield_tmp_vec2_(3 * size_phi_, 0., 0)
  , phi_tmp_vec_(size_phi_, 0., 0)
  , phi_tmp_vec2_(size_phi_, 0., 0)
  , u_tmp_vec_(size_u_, 0., 0)
  , P_tmp_vec_(size_P_, 0., 0)
  , u_tmp_(u_space_)
  , u_tmp_local_(std::make_shared<PerThreadVectorLocalFunc>())
  , P_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
  , Pnat_tmp_local_(std::make_shared<PerThreadVectorLocalFuncs>(num_cells_))
  , phi_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
  , phinat_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
  , mu_tmp_local_(std::make_shared<PerThreadScalarLocalFuncs>(num_cells_))
{
  std::cout << dx_ << ", " << epsilon_ << std::endl;

  /************************** create and project initial values*****************************************
   ************************** we only need initial values for P and phi *******************************/

  std::shared_ptr<const XT::Functions::FunctionInterface<d, d>> u_initial_func;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> phi_initial_funcs;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d, d>>> P_initial_funcs;

  // interpolate initial and boundary values
  if (testcase == "single_cell" || testcase == "channel") {
    // Initially, cell is circular with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 30]^2
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(15, 15) - x|.
    FieldVector<double, d> center{upper_right_[0] / 2., upper_right_[1] / 2.};
    auto r = [center](const auto& xr) { return 5.0 - (center - xr).two_norm(); };
    phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
        50,
        /*evaluate=*/
        [r, epsilon = epsilon_](const auto& x, const auto& /*param*/) {
          return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
        },
        /*name=*/"phi_initial"));

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
  } else if (testcase == "two_cells") {
    FieldVector<double, d> center1{15, 15};
    FieldVector<double, d> center2{35, 35};
    auto r1 = [center1](const auto& xr) { return 4.0 - (center1 - xr).two_norm(); };
    auto r2 = [center2](const auto& xr) { return 4.0 - (center2 - xr).two_norm(); };
    phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
        50,
        /*evaluate=*/
        [r = r1, epsilon = epsilon_](const auto& x, const auto& /*param*/) {
          return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
        },
        /*name=*/"phi1_initial"));
    phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
        50,
        /*evaluate=*/
        [r = r2, epsilon = epsilon_](const auto& x, const auto& /*param*/) {
          return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
        },
        /*name=*/"phi2_initial"));

    // initial condition for P is (1,0) + \delta where \delta(x) is vector-valued with random entries following an
    // uniform distribution on the interval [-0.05, 0.05]; restrict to cytoplasm by multiplying with (\phi + 1)/2
    std::srand(1); // set seed for std::rand to 1
    P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
        50,
        /*evaluate=*/
        [& phi1_initial = phi_initial_funcs[0]](const auto& x, const auto& param) {
          // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto ret = FieldVector<double, d>({1. + rand1, 0. +
          // rand2});
          auto ret = FieldVector<double, d>({1., 0.});
          ret *= (phi1_initial->evaluate(x, param) + 1.) / 2.;
          return ret;
        },
        /*name=*/"P1_initial"));
    P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
        50,
        /*evaluate=*/
        [& phi2_initial = phi_initial_funcs[1]](const auto& x, const auto& param) {
          // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto ret = FieldVector<double, d>({1. + rand1, 0. +
          // rand2});
          auto ret = FieldVector<double, d>({1., 0.});
          ret *= (phi2_initial->evaluate(x, param) + 1.) / 2.;
          return ret;
        },
        /*name=*/"P2_initial"));

    u_initial_func = std::make_shared<const XT::Functions::ConstantFunction<d, d>>(0.);
  } else {
    DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  }

  /*************************************************************************************************
   ******************************* create variables, set initial values ****************************
   *************************************************************************************************/

  // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0, Neumann boundary
  // conditions for the other variables
  XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
  default_interpolation(*u_initial_func, u_);
  // create system and temporary vectors, DiscreteFunctions, etc.
  for (size_t kk = 0; kk < num_cells_; kk++) {
    P_view_.emplace_back(ofield_vectors_[kk], 0, size_P_);
    Pnat_view_.emplace_back(ofield_vectors_[kk], size_P_, 2 * size_P_);
    phi_view_.emplace_back(pfield_vectors_[kk], 0, size_phi_);
    phinat_view_.emplace_back(pfield_vectors_[kk], size_phi_, 2 * size_phi_);
    mu_view_.emplace_back(pfield_vectors_[kk], 2 * size_phi_, 3 * size_phi_);
    const auto kk_str = XT::Common::to_string(kk);
    P_.emplace_back(make_discrete_function(P_space_, P_view_[kk], "P_" + kk_str));
    Pnat_.emplace_back(make_discrete_function(P_space_, Pnat_view_[kk], "Pnat_" + kk_str));
    phi_.emplace_back(make_discrete_function(phi_space_, phi_view_[kk], "phi_" + kk_str));
    phinat_.emplace_back(make_discrete_function(phi_space_, phinat_view_[kk], "phinat_" + kk_str));
    mu_.emplace_back(make_discrete_function(phi_space_, mu_view_[kk], "mu_" + kk_str));
    P_tmp_.emplace_back(P_space_);
    Pnat_tmp_.emplace_back(P_space_);
    phi_tmp_.emplace_back(phi_space_);
    phinat_tmp_.emplace_back(phi_space_);
    mu_tmp_.emplace_back(phi_space_);
    default_interpolation(*phi_initial_funcs[kk], phi_[kk]);
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
  p_basis_integrated_functional.append(
      LocalElementIntegralFunctional<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>().with_ansatz(one_function)));
  B_stokes_op.append(p_basis_integrated_functional);

  // Dirichlet constrainst for u
  A_stokes_op_->append(u_dirichlet_constraints_);
  // assemble everything
  A_stokes_op_->append(B_stokes_op);
  A_stokes_op_->append(M_p_stokes_op);
  std::cout << "Assembling stokes..." << std::flush;
  A_stokes_op_->assemble(use_tbb_);
  std::cout << "done" << std::endl;

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

  std::cout << "Setting up stokes solver..." << std::flush;
  S_colmajor_ = S_stokes_.backend();
  S_colmajor_.makeCompressed();
  stokes_solver_->analyzePattern(S_colmajor_);
  stokes_solver_->factorize(S_colmajor_);
  std::cout << "done" << std::endl;

  /*************************************************************************************************
   ************************************ Orientationfield *******************************************
   *************************************************************************************************/
  // calculate M_{ij} as \int \psi_i phi_j
  M_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(1.)));
  C_ofield_elliptic_part_ *= 0.;
  MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, P_space_, P_space_, C_ofield_elliptic_part_);
  ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
  M_ofield_op_->append(ofield_elliptic_op);

  std::cout << "Assembling ofield..." << std::flush;
  M_ofield_op_->assemble(use_tbb_);
  std::cout << "done" << std::endl;
  std::cout << "Setting up ofield_solver..." << std::flush;
  ofield_solver_.setup();
  std::cout << "done" << std::endl;

  /*************************************************************************************************
   **************************************** Phasefield *********************************************
   *************************************************************************************************/

  MatrixOperator<MatrixType, PGV, 1> M_pfield_op(grid_view_, phi_space_, phi_space_, M_pfield_);
  M_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(1.)));
  MatrixOperator<MatrixType, PGV, 1> M_ell_pfield_op(grid_view_, phi_space_, phi_space_, M_ell_pfield_);
  M_ell_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
  M_pfield_op.append(M_ell_pfield_op);
  std::cout << "Assembling pfield..." << std::flush;
  M_pfield_op.assemble(use_tbb_);
  std::cout << "done" << std::endl;
  std::cout << "Setting up pfield_solver..." << std::flush;
  pfield_solver_.setup();
  std::cout << "done" << std::endl;
} // constructor

size_t CellModelSolver::num_cells() const
{
  return num_cells_;
}

bool CellModelSolver::linear() const
{
  return false;
}

bool CellModelSolver::finished() const
{
  return XT::Common::FloatCmp::eq(t_end_, t_);
}

//******************************************************************************************************************
//********************************* Solve methods for whole system of equations ************************************
//******************************************************************************************************************

/**
 * Solves whole system of equations using the values stored in stokes_vector_, ofield_vectors_ and pfield_vectors_ as
 *initial values. Returns the whole trajectory, i.e., ret[i] contains the results in the i-th timestep. The first
 *num_cells_ entries of ret[i] correspond to the phasefield for each cell, the next num_cells entries are the
 *orientation field vectors, and the last one is the stokes vector. dt: Time step length. write: Whether to write .vtu
 *and .txt files. write_step: Time interval at which results should be written. If negative, all steps are  written.
 *Ignored if write = false. filename: Prefix for .vtu and .txt files. Ignored if write = false. subsampling: Whether to
 *use subsampling for visualization. Ignored if write = false.
 **/
std::vector<std::vector<CellModelSolver::VectorType>>
CellModelSolver::solve(const bool write, const double write_step, const std::string filename, const bool subsampling)
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
  assert(XT::Common::FloatCmp::ge(t_end_, t_));
  double next_save_time = t_ + write_step > t_end_ ? t_end_ : t_ + write_step;
  size_t save_step_counter = 1;

  while (XT::Common::FloatCmp::lt(t_, t_end_)) {
    double max_dt = dt_;
    // match saving times and t_end_ exactly
    if (XT::Common::FloatCmp::gt(t_ + dt_, t_end_))
      max_dt = t_end_ - t_;
    double actual_dt = std::min(dt_, max_dt);

    // do a timestep
    std::cout << "Current time: " << t_ << std::endl;
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      prepare_pfield_operator(kk);
      ret[kk].push_back(apply_inverse_pfield_operator(ret[kk].back(), kk));
      set_pfield_vec(kk, ret[kk].back());
      std::cout << "Pfield " << kk << " done" << std::endl;
      prepare_ofield_operator(kk);
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
      if (write_step < 0. || XT::Common::FloatCmp::ge(t_, next_save_time)) {
        visualize(filename, save_step_counter, t_, subsampling);
        next_save_time += write_step;
        ++save_step_counter;
      }
    }
  } // while (t_ < t_end_)
  return ret;
}

// Like solve, but only computes and returns the next n timesteps
std::vector<std::vector<CellModelSolver::VectorType>> CellModelSolver::next_n_timesteps(const size_t n)
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

  assert(XT::Common::FloatCmp::ge(t_end_, t_));

  // implicit Euler timestepping
  while (XT::Common::FloatCmp::lt(t_, t_end_) && count < n) {
    double max_dt = dt_;
    // match saving times and t_end_ exactly
    if (XT::Common::FloatCmp::gt(t_ + dt_, t_end_))
      max_dt = t_end_ - t_;
    double actual_dt = std::min(dt_, max_dt);

    // do a timestep
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      prepare_pfield_operator(kk);
      pfield_vectors_[kk] = apply_inverse_pfield_operator(pfield_vectors_[kk], kk);
      ret[kk].push_back(pfield_vectors_[kk]);
      // std::cout << "Pfield " << kk << " done" << std::endl;
      prepare_ofield_operator(kk);
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
  VectorType ret(2 * size_P_);
  ConstVectorViewType P_view(vec, 0, size_P_);
  ConstVectorViewType Pnat_view(vec, size_P_, 2 * size_P_);
  VectorViewType P_ret_view(ret, 0, size_P_);
  VectorViewType Pnat_ret_view(ret, size_P_, 2 * size_P_);
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
  auto vtk_writer = XT::Functions::internal::create_vtkwriter(phi_space_.grid_view(), subsampling);
  VectorType phi_vec(size_phi_);
  for (size_t ii = 0; ii < size_phi_; ++ii)
    phi_vec[ii] = vec[ii];
  const ConstVectorViewType phinat_vec(vec, size_phi_, 2 * size_phi_);
  const ConstVectorViewType mu_vec(vec, 2 * size_phi_, 3 * size_phi_);
  const auto phi_func = make_discrete_function(phi_space_, phi_vec, "phi");
  const auto phinat_func = make_discrete_function(phi_space_, phinat_vec, "phinat");
  const auto mu_func = make_discrete_function(phi_space_, mu_vec, "mu");
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phi_func);
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phinat_func);
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, mu_func);
  XT::Functions::internal::write_visualization(*vtk_writer, filename);
} // void visualize_pfield(...)

// Visualizes given vector as orientation field finite element vector
void CellModelSolver::visualize_ofield(const std::string& filename, const VectorType& vec, const bool subsampling) const
{
  auto vtk_writer = XT::Functions::internal::create_vtkwriter(P_space_.grid_view(), subsampling);
  const ConstVectorViewType P_vec(vec, 0, size_P_);
  const ConstVectorViewType Pnat_vec(vec, size_P_, 2 * size_P_);
  const auto P_func = make_discrete_function(P_space_, P_vec, "P");
  const auto Pnat_func = make_discrete_function(P_space_, Pnat_vec, "Pnat");
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, P_func);
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, Pnat_func);
  XT::Functions::internal::write_visualization(*vtk_writer, filename);
} // void visualize_ofield(...)

// Visualizes given vector as stokes finite element vector
void CellModelSolver::visualize_stokes(const std::string& filename, const VectorType& vec, const bool subsampling) const
{
  auto vtk_writer = XT::Functions::internal::create_vtkwriter(u_.space().grid_view(), subsampling);
  const ConstVectorViewType u_vec(vec, 0, size_u_);
  const ConstVectorViewType p_vec(vec, size_u_, size_u_ + size_p_);
  const auto u_func = make_discrete_function(u_space_, u_vec, "u");
  const auto p_func = make_discrete_function(p_space_, p_vec, "p");
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, u_func);
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, p_func);
  XT::Functions::internal::write_visualization(*vtk_writer, filename);
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
  auto vtk_writer = XT::Functions::internal::create_vtkwriter(u_.space().grid_view(), subsampling);
  std::string postfix = "_" + XT::Common::to_string(step);
  if (vtu) {
    XT::Functions::internal::add_to_vtkwriter(*vtk_writer, u_);
    XT::Functions::internal::add_to_vtkwriter(*vtk_writer, p_);
    std::vector<VectorType> phi_vec(num_cells_, VectorType(size_phi_, 0.));
    std::vector<std::shared_ptr<ConstDiscreteFunctionType>> phi_funcs(num_cells_);
    for (size_t kk = 0; kk < num_cells_; ++kk) {
      // for (size_t ii = 0; ii < size_phi_; ++ii)
      // phi_vec[kk][ii] = phi_[kk].dofs().vector()[ii];
      // phi_funcs[kk] =
      // std::make_shared<ConstDiscreteFunctionType>(phi_space_, phi_vec[kk], "phi_" + XT::Common::to_string(kk));
      XT::Functions::internal::add_to_vtkwriter(*vtk_writer, P_[kk]);
      XT::Functions::internal::add_to_vtkwriter(*vtk_writer, Pnat_[kk]);
      // XT::Functions::internal::add_to_vtkwriter(*vtk_writer, *phi_funcs[kk]);
      XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phi_[kk]);
      XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phinat_[kk]);
      XT::Functions::internal::add_to_vtkwriter(*vtk_writer, mu_[kk]);
      // XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, *phi_funcs[kk]);
      XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, phi_[kk]);
      XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, phinat_[kk]);
      XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, mu_[kk]);
    }
    XT::Functions::internal::write_visualization(*vtk_writer, prefix + postfix);
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

// Sets given dofs of stokes vector to values
void CellModelSolver::set_stokes_vec_dofs(const std::vector<R>& values, const std::vector<size_t>& dofs)
{
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    stokes_vector_.set_entry(dofs[ii], values[ii]);
}

// Sets orientation field vector belonging to cell to pfield_vec
void CellModelSolver::set_ofield_vec(const size_t cell, const VectorType& ofield_vec)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(ofield_vec.size() != 2 * size_P_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
  ofield_vectors_[cell] = ofield_vec;
}

// Sets given dofs of orientation field vector belonging to cell to values
void CellModelSolver::set_ofield_vec_dofs(const size_t cell,
                                          const std::vector<R>& values,
                                          const std::vector<size_t>& dofs)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    ofield_vectors_[cell].set_entry(dofs[ii], values[ii]);
}

// Sets phasefield vector belonging to cell to pfield_vec
void CellModelSolver::set_pfield_vec(const size_t cell, const VectorType& pfield_vec)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(pfield_vec.size() != 3 * size_phi_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
  pfield_vectors_[cell] = pfield_vec;
}

// Sets given dofs of phasefield vector belonging to cell to values
void CellModelSolver::set_pfield_vec_dofs(const size_t cell,
                                          const std::vector<R>& values,
                                          const std::vector<size_t>& dofs)
{
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    pfield_vectors_[cell].set_entry(dofs[ii], values[ii]);
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

void CellModelSolver::prepare_stokes_operator(const bool restricted)
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  for (size_t kk = 0; kk < num_cells_; kk++) {
    phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
    phinat_tmp_[kk].dofs().vector() = phinat_[kk].dofs().vector();
    P_tmp_[kk].dofs().vector() = P_[kk].dofs().vector();
    Pnat_tmp_[kk].dofs().vector() = Pnat_[kk].dofs().vector();
  }
  assemble_stokes_rhs(restricted);
}

void CellModelSolver::prepare_ofield_operator(const size_t cell, const bool restricted)
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
  phi_tmp_[cell].dofs().vector() = phi_[cell].dofs().vector();
  ofield_jac_linear_op_.prepare(cell, restricted);
  assemble_ofield_rhs(cell, restricted);
  assemble_ofield_linear_jacobian(cell, restricted);
  if (!restricted)
    ofield_solver_.prepare(cell, restricted);
}

void CellModelSolver::prepare_pfield_operator(const size_t cell, const bool restricted)
{
  u_tmp_.dofs().vector() = u_.dofs().vector();
  P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
  for (size_t kk = 0; kk < num_cells_; kk++) {
    phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
    mu_tmp_[kk].dofs().vector() = mu_[kk].dofs().vector();
  }
  assemble_pfield_rhs(cell, restricted);
  assemble_pfield_linear_jacobian(cell, restricted);
  pfield_jac_linear_op_.prepare(cell, restricted);
  if (!restricted)
    pfield_solver_.prepare(cell, restricted);
}

void CellModelSolver::compute_restricted_stokes_dofs(const std::vector<size_t>& range_dofs)
{
  if (!stokes_deim_range_dofs_ || *stokes_deim_range_dofs_ != range_dofs) {
    // We need to keep the original range_dofs which is unordered and may contain duplicates, as the restricted
    // operator will return exactly these dofs. For computations, however, we often need unique dofs.
    stokes_deim_range_dofs_ = std::make_shared<std::vector<size_t>>(range_dofs);
    auto& unique_range_dofs = stokes_deim_unique_range_dofs_;
    get_unique_deim_dofs(unique_range_dofs, range_dofs, size_u_ + size_p_);
    // sort output into dofs belonging to u and p
    auto& u_range_dofs = u_deim_range_dofs_;
    auto& p_range_dofs = p_deim_range_dofs_;
    u_range_dofs.clear();
    p_range_dofs.clear();
    for (const auto& dof : unique_range_dofs) {
      if (dof < size_u_)
        u_range_dofs.push_back(dof);
      else
        p_range_dofs.push_back(dof);
    }
    subtract_from_dofs(p_range_dofs, size_u_);
    // store all entities that contain an range dof
    const auto& u_mapper = u_space_.mapper();
    const auto& p_mapper = p_space_.mapper();
    DynamicVector<size_t> global_indices;
    stokes_deim_entities_.clear();
    for (const auto& entity : Dune::elements(grid_view_))
      maybe_add_entity_stokes(
          entity, global_indices, u_range_dofs, p_range_dofs, stokes_deim_entities_, u_mapper, p_mapper);
    // get source dofs corresponding to range dofs
    auto& source_dofs = stokes_deim_source_dofs_;
    source_dofs.clear();
    source_dofs.resize(3);
    get_deim_source_dofs(source_dofs[2], create_stokes_pattern(u_space_, p_space_), unique_range_dofs);
    // add phi, phinat, P and Pnat source dofs for all input entities
    const auto& phi_mapper = phi_space_.mapper();
    const auto& P_mapper = P_space_.mapper();
    for (const auto& entity : stokes_deim_entities_) {
      global_indices.resize(phi_mapper.local_size(entity));
      phi_mapper.global_indices(entity, global_indices);
      source_dofs[0].insert(source_dofs[0].end(), global_indices.begin(), global_indices.end());
      for (auto& index : global_indices)
        index += size_phi_;
      source_dofs[0].insert(source_dofs[0].end(), global_indices.begin(), global_indices.end());
      global_indices.resize(P_mapper.local_size(entity));
      P_mapper.global_indices(entity, global_indices);
      source_dofs[1].insert(source_dofs[1].end(), global_indices.begin(), global_indices.end());
      for (auto& index : global_indices)
        index += size_P_;
      source_dofs[1].insert(source_dofs[1].end(), global_indices.begin(), global_indices.end());
    } // entities
    // sort and remove duplicate entries
    sort_and_remove_duplicates_in_deim_source_dofs(source_dofs);
  } // if (not already computed)
} // void compute_restricted_stokes_dofs(...)

void CellModelSolver::compute_restricted_ofield_dofs(const std::vector<size_t>& range_dofs, const size_t cell)
{
  if (!ofield_deim_range_dofs_[cell] || *ofield_deim_range_dofs_[cell] != range_dofs) {
    // We need to keep the original range_dofs which is unordered and may contain duplicates, as the restricted
    // operator will return exactly these dofs. For computations, however, we often need unique dofs.
    ofield_deim_range_dofs_[cell] = std::make_shared<std::vector<size_t>>(range_dofs);
    auto& unique_range_dofs = ofield_deim_unique_range_dofs_[cell];
    get_unique_deim_dofs(unique_range_dofs, range_dofs, 2 * size_P_);
    // sort output into dofs belonging to P and Pnat
    auto& P_range_dofs = P_deim_range_dofs_[cell];
    auto& Pnat_range_dofs = Pnat_deim_range_dofs_[cell];
    P_range_dofs.clear();
    Pnat_range_dofs.clear();
    for (const auto& dof : unique_range_dofs) {
      if (dof < size_P_)
        P_range_dofs.push_back(dof);
      else
        Pnat_range_dofs.push_back(dof);
    }
    subtract_from_dofs(Pnat_range_dofs, size_P_);
    // store all entities that contain a range dof
    const auto& P_mapper = P_space_.mapper();
    get_deim_entities(P_mapper, ofield_deim_entities_[cell], unique_range_dofs, size_P_);

    // get source dofs corresponding to range dofs from system matrix pattern
    auto& source_dofs = ofield_deim_source_dofs_[cell];
    source_dofs.clear();
    source_dofs.resize(3);
    get_deim_source_dofs(
        source_dofs[1], ofield_solver_.system_matrix_pattern(ofield_submatrix_pattern_), unique_range_dofs);

    // add phi and u source dofs for all input entities
    DynamicVector<size_t> global_indices;
    const auto& phi_mapper = phi_space_.mapper();
    const auto& u_mapper = u_space_.mapper();
    for (const auto& entity : ofield_deim_entities_[cell]) {
      global_indices.resize(phi_mapper.local_size(entity));
      phi_mapper.global_indices(entity, global_indices);
      source_dofs[0].insert(source_dofs[0].end(), global_indices.begin(), global_indices.end());
      global_indices.resize(u_mapper.local_size(entity));
      u_mapper.global_indices(entity, global_indices);
      source_dofs[2].insert(source_dofs[2].end(), global_indices.begin(), global_indices.end());
    } // entities

    // sort and remove duplicate entries
    sort_and_remove_duplicates_in_deim_source_dofs(source_dofs);

    Pnat_deim_source_dofs_begin_[cell] =
        std::lower_bound(source_dofs[1].begin(), source_dofs[1].end(), size_P_) - source_dofs[1].begin();
  } // if (not already computed)
} // void compute_restricted_ofield_dofs(...)

void CellModelSolver::compute_restricted_pfield_dofs(const std::vector<size_t>& range_dofs, const size_t cell)
{
  if (!pfield_deim_range_dofs_[cell] || *pfield_deim_range_dofs_[cell] != range_dofs) {
    // We need to keep the original range_dofs which is unordered and may contain duplicates, as the restricted
    // operator will return exactly these dofs. For computations, however, we often need unique dofs.
    pfield_deim_range_dofs_[cell] = std::make_shared<std::vector<size_t>>(range_dofs);
    auto& unique_range_dofs = pfield_deim_unique_range_dofs_[cell];
    get_unique_deim_dofs(unique_range_dofs, range_dofs, 3 * size_phi_);
    // sort output into dofs belonging to phi, phinat and mu
    auto& phi_range_dofs = phi_deim_range_dofs_[cell];
    auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
    auto& mu_range_dofs = mu_deim_range_dofs_[cell];
    auto& phinat_mu_range_dofs = both_mu_and_phi_deim_range_dofs_[cell];
    phi_range_dofs.clear();
    phinat_range_dofs.clear();
    mu_range_dofs.clear();
    for (const auto& dof : unique_range_dofs) {
      if (dof < size_phi_)
        phi_range_dofs.push_back(dof);
      else if (dof < 2 * size_phi_)
        phinat_range_dofs.push_back(dof);
      else
        mu_range_dofs.push_back(dof);
    }
    subtract_from_dofs(phinat_range_dofs, size_phi_);
    subtract_from_dofs(mu_range_dofs, 2 * size_phi_);
    // get dofs that are both in phinat and mu dofs
    phinat_mu_range_dofs.clear();
    for (const auto& dof : phinat_range_dofs)
      if (std::find(mu_range_dofs.begin(), mu_range_dofs.end(), dof) != mu_range_dofs.end())
        phinat_mu_range_dofs.push_back(dof);
    // store all entities that contain a range dof
    get_deim_entities(phi_space_.mapper(), pfield_deim_entities_[cell], unique_range_dofs, size_phi_);

    // get source dofs corresponding to range dofs
    auto& source_dofs = pfield_deim_source_dofs_[cell];
    source_dofs.clear();
    source_dofs.resize(3);
    get_deim_source_dofs(
        source_dofs[0], pfield_solver_.system_matrix_pattern(pfield_submatrix_pattern_), unique_range_dofs);

    // add P and u source dofs for all input entities
    DynamicVector<size_t> global_indices;
    const auto& P_mapper = P_space_.mapper();
    const auto& u_mapper = u_space_.mapper();
    for (const auto& entity : pfield_deim_entities_[cell]) {
      global_indices.resize(P_mapper.local_size(entity));
      P_mapper.global_indices(entity, global_indices);
      source_dofs[1].insert(source_dofs[1].end(), global_indices.begin(), global_indices.end());
      global_indices.resize(u_mapper.local_size(entity));
      u_mapper.global_indices(entity, global_indices);
      source_dofs[2].insert(source_dofs[2].end(), global_indices.begin(), global_indices.end());
    } // entities

    // sort and remove duplicate entries
    sort_and_remove_duplicates_in_deim_source_dofs(source_dofs);

    phinat_deim_source_dofs_begin_[cell] =
        std::lower_bound(source_dofs[0].begin(), source_dofs[0].end(), size_phi_) - source_dofs[0].begin();
    mu_deim_source_dofs_begin_[cell] =
        std::lower_bound(source_dofs[0].begin(), source_dofs[0].end(), 2 * size_phi_) - source_dofs[0].begin();
  } // if (not already computed)
} // void compute_restricted_pfield_dofs(...)

// appends entity to input_entities if one of its global_indices is in range_dofs
void CellModelSolver::maybe_add_entity(const E& entity,
                                       DynamicVector<size_t>& global_indices,
                                       const std::vector<size_t>& range_dofs,
                                       std::vector<E>& input_entities,
                                       const MapperInterface<PGV>& mapper,
                                       const size_t subvector_size) const
{
  global_indices.resize(mapper.local_size(entity));
  mapper.global_indices(entity, global_indices);
  for (const auto& output_dof : range_dofs) {
    const size_t dof = output_dof % subvector_size;
    for (size_t jj = 0; jj < global_indices.size(); ++jj) {
      if (global_indices[jj] == dof) {
        input_entities.push_back(entity);
        return;
      }
    } // jj
  } // dof
}

// appends entity to input_entities if one of its global_indices is in range_dofs
void CellModelSolver::maybe_add_entity_stokes(const E& entity,
                                              DynamicVector<size_t>& global_indices,
                                              const std::vector<size_t>& u_range_dofs,
                                              const std::vector<size_t>& p_range_dofs,
                                              std::vector<E>& stokes_deim_entities,
                                              const MapperInterface<PGV>& u_mapper,
                                              const MapperInterface<PGV>& p_mapper)
{
  global_indices.resize(u_mapper.local_size(entity));
  u_mapper.global_indices(entity, global_indices);
  for (const auto& dof : u_range_dofs) {
    for (size_t jj = 0; jj < global_indices.size(); ++jj) {
      if (global_indices[jj] == dof) {
        stokes_deim_entities.push_back(entity);
        return;
      }
    } // jj
  } // dof
  global_indices.resize(p_mapper.local_size(entity));
  p_mapper.global_indices(entity, global_indices);
  for (const auto& dof : p_range_dofs) {
    for (size_t jj = 0; jj < global_indices.size(); ++jj) {
      if (global_indices[jj] == dof) {
        stokes_deim_entities.push_back(entity);
        return;
      }
    } // jj
  } // dof
}

void CellModelSolver::subtract_from_dofs(std::vector<size_t>& dofs, const size_t size) const
{
  for (auto& dof : dofs) {
    assert(dof >= size);
    dof -= size;
  }
}

void CellModelSolver::get_deim_entities(const MapperInterface<PGV>& mapper,
                                        std::vector<E>& deim_entities,
                                        const std::vector<size_t>& unique_range_dofs,
                                        const size_t subvector_size) const
{
  DynamicVector<size_t> global_indices;
  deim_entities.clear();
  for (const auto& entity : Dune::elements(grid_view_)) {
    maybe_add_entity(entity, global_indices, unique_range_dofs, deim_entities, mapper, subvector_size);
  } // entities
}

void CellModelSolver::get_unique_deim_dofs(std::vector<size_t>& unique_range_dofs,
                                           const std::vector<size_t>& range_dofs,
                                           DXTC_DEBUG_ONLY const size_t max_dof_value) const
{
  unique_range_dofs = range_dofs;
  std::sort(unique_range_dofs.begin(), unique_range_dofs.end());
  unique_range_dofs.erase(std::unique(unique_range_dofs.begin(), unique_range_dofs.end()), unique_range_dofs.end());
  // check that all dofs are valid
  DEBUG_THROW_IF(!XT::Common::transform_reduce(
                     unique_range_dofs.begin(),
                     unique_range_dofs.end(),
                     true,
                     [](const bool& a, const bool& b) { return a && b; },
                     [size = max_dof_value](const size_t& a) { return a < size; }),
                 XT::Common::Exceptions::wrong_input_given,
                 "At least one of the given output DoFs is too large!");
}

void CellModelSolver::get_deim_source_dofs(std::vector<size_t>& source_dofs,
                                           const XT::LA::SparsityPatternDefault& pattern,
                                           const std::vector<size_t>& unique_range_dofs) const
{
  for (const auto& dof : unique_range_dofs) {
    const auto& new_source_dofs = pattern.inner(dof);
    source_dofs.insert(source_dofs.end(), new_source_dofs.begin(), new_source_dofs.end());
  }
}

void CellModelSolver::sort_and_remove_duplicates_in_deim_source_dofs(
    std::vector<std::vector<size_t>>& source_dofs) const
{
  // check that source_dofs has entries for pfield, ofield and stokes
  assert(source_dofs.size() == 3);
  std::vector<size_t> sizes{3 * size_phi_, 2 * size_P_, size_u_ + size_p_};
  for (size_t ii = 0; ii < 3; ++ii) {
    std::sort(source_dofs[ii].begin(), source_dofs[ii].end());
    source_dofs[ii].erase(std::unique(source_dofs[ii].begin(), source_dofs[ii].end()), source_dofs[ii].end());
    // check that dofs are valid
    DEBUG_THROW_IF(!XT::Common::transform_reduce(
                       source_dofs[ii].begin(),
                       source_dofs[ii].end(),
                       true,
                       [](const bool& a, const bool& b) { return a && b; },
                       [size = sizes[ii]](const size_t& a) { return a < size; }),
                   XT::Common::Exceptions::wrong_input_given,
                   "At least one of the given output DoFs is too large!");
  }
}

//******************************************************************************************************************
//*********************************************** Apply operators **************************************************
//******************************************************************************************************************

CellModelSolver::VectorType CellModelSolver::apply_stokes_operator(const VectorType& y, const bool restricted)
{
  return apply_stokes_helper(y, restricted, false);
}

// Applies stokes operator (applies the F if Stokes equation is F(y) = 0)
CellModelSolver::VectorType
CellModelSolver::apply_stokes_helper(const VectorType& y, const bool restricted, const bool jacobian)
{
  const auto& unique_range_dofs = stokes_deim_unique_range_dofs_;
  auto& source = stokes_tmp_vec_;
  auto& residual = stokes_tmp_vec2_;
  // copy values to high-dimensional vector
  if (restricted)
    copy_ld_to_hd_vec(stokes_deim_source_dofs_[2], y, source);
  else
    source = y;
  const auto mv = mv_func<VectorType, VectorType>(restricted);
  mv(S_stokes_, source, residual, unique_range_dofs);
  if (!jacobian) {
    const auto sub = sub_func<VectorType, EigenVectorType>(restricted);
    sub(residual, stokes_rhs_vector_, unique_range_dofs);
  }
  if (restricted) {
    const auto& range_dofs = *stokes_deim_range_dofs_;
    VectorType ret(range_dofs.size());
    for (size_t ii = 0; ii < range_dofs.size(); ++ii)
      ret[ii] = residual[range_dofs[ii]];
    return ret;
  } else {
    return residual;
  }
}

void CellModelSolver::assemble_nonlinear_part_of_ofield_residual(VectorType& residual,
                                                                 const size_t cell,
                                                                 const bool restricted)
{
  // nonlinear part
  VectorViewType res1_vec(residual, size_P_, 2 * size_P_);
  auto nonlinear_res_functional = make_vector_functional(P_space_, res1_vec);
  XT::Functions::GenericGridFunction<E, d, 1> nonlinear_res_pf(
      /*order = */ 3 * P_space_.max_polorder(),
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
      LocalElementProductScalarWeightIntegrand<E, d>().with_ansatz(nonlinear_res_pf)));
  if (!restricted)
    nonlinear_res_functional.assemble(use_tbb_);
  else
    nonlinear_res_functional.assemble_range(ofield_deim_entities_[cell]);
}

// Applies cell-th orientation field operator (applies F if the orientation field equation is F(y) = 0)
CellModelSolver::VectorType
CellModelSolver::apply_ofield_operator(const VectorType& y, const size_t cell, const bool restricted)
{
  return apply_ofield_helper(y, cell, restricted, false);
}

CellModelSolver::VectorType
CellModelSolver::apply_ofield_helper(const VectorType& y, const size_t cell, const bool restricted, const bool jacobian)
{
  auto& source = ofield_tmp_vec_;
  auto& residual = ofield_tmp_vec2_;
  // copy values to high-dimensional vector
  if (restricted)
    copy_ld_to_hd_vec(ofield_deim_source_dofs_[cell][1], y, source);
  else
    source = y;
  // linear part
  ofield_jac_linear_op_.apply(source, residual);
  if (!jacobian) {
    // subtract rhs
    const auto& unique_range_dofs = ofield_deim_unique_range_dofs_[cell];
    const auto sub = sub_func<VectorType>(restricted);
    sub(residual, ofield_rhs_vector_, unique_range_dofs);
  }
  // nonlinear part
  if (jacobian) {
    const auto& Pnat_range_dofs = Pnat_deim_range_dofs_[cell];
    VectorViewType range_Pnat(residual, size_P_, 2 * size_P_);
    const ConstVectorViewType source_P(source, 0, size_P_);
    const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
    const auto add = add_func<VectorViewType, VectorType>(restricted);
    mv(C_ofield_nonlinear_part_, source_P, P_tmp_vec_, Pnat_range_dofs);
    add(range_Pnat, P_tmp_vec_, Pnat_range_dofs);
  } else {
    fill_tmp_ofield(cell, source, restricted);
    assemble_nonlinear_part_of_ofield_residual(residual, cell, restricted);
  }
  if (restricted) {
    const auto& range_dofs = *ofield_deim_range_dofs_[cell];
    VectorType ret(range_dofs.size());
    for (size_t ii = 0; ii < range_dofs.size(); ++ii)
      ret[ii] = residual[range_dofs[ii]];
    return ret;
  } else {
    return residual;
  }
}

void CellModelSolver::update_ofield_parameters(const double Pa, const size_t cell, const bool restricted)
{
  DUNE_THROW_IF(
      cell != 0, Dune::NotImplemented, "This may not work for several cells, check before removing this line!");
  // Pa may have been set to a new value already (via update_pfield_parameters)
  if (XT::Common::FloatCmp::ne(Pa, last_ofield_Pa_)) {
    std::cout << "Ofield params updated, old Pa = " << last_ofield_Pa_ << ", new Pa = " << Pa << std::endl;
    Pa_ = Pa;
    last_ofield_Pa_ = Pa_;
    const auto& P_range_dofs = P_deim_range_dofs_[cell];
    set_mat_to_zero(C_ofield_elliptic_part_, restricted, ofield_submatrix_pattern_, P_range_dofs);
    MatrixOperator<MatrixType, PGV, d> ofield_elliptic_op(grid_view_, P_space_, P_space_, C_ofield_elliptic_part_);
    ofield_elliptic_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
    if (!restricted)
      ofield_elliptic_op.assemble(use_tbb_);
    else
      ofield_elliptic_op.assemble_range(ofield_deim_entities_[cell]);
  }
}

// Applies cell-th phase field operator (applies F if phase field equation is F(y) = 0)
CellModelSolver::VectorType
CellModelSolver::apply_pfield_operator(const VectorType& y, const size_t cell, const bool restricted)
{
  return apply_pfield_helper(y, cell, restricted, false);
}

CellModelSolver::VectorType
CellModelSolver::apply_pfield_helper(const VectorType& y, const size_t cell, const bool restricted, const bool jacobian)
{
  auto& source = pfield_tmp_vec_;
  auto& residual = pfield_tmp_vec2_;
  // copy values to high-dimensional vector
  if (restricted)
    copy_ld_to_hd_vec(pfield_deim_source_dofs_[cell][0], y, source);
  else
    source = y;
  // linear part
  pfield_jac_linear_op_.apply(source, residual);
  if (!jacobian) {
    // subtract rhs
    const auto& unique_range_dofs = pfield_deim_unique_range_dofs_[cell];
    const auto sub = sub_func<VectorType>(restricted);
    sub(residual, pfield_rhs_vector_, unique_range_dofs);
  }
  // nonlinear part
  if (jacobian) {
    const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
    const auto& mu_range_dofs = mu_deim_range_dofs_[cell];
    VectorViewType range_phinat(residual, size_phi_, 2 * size_phi_);
    VectorViewType range_mu(residual, 2 * size_phi_, 3 * size_phi_);
    const ConstVectorViewType source_phi(source, 0, size_phi_);
    const ConstVectorViewType source_mu(source, 2 * size_phi_, 3 * size_phi_);
    // nonlinear_part
    auto& tmp_vec = phi_tmp_vec_;
    const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
    const auto vector_axpy = vector_axpy_func<VectorViewType, VectorType>(restricted);
    const auto scal = scal_func<VectorType>(restricted);
    const auto add = add_func<VectorViewType, VectorType>(restricted);
    // apply missing parts of J (including the linear 1./Ca_ part)
    mv(M_nonlin_pfield_, source_mu, tmp_vec, phinat_range_dofs);
    vector_axpy(range_phinat, 1. / (Be_ * std::pow(epsilon_, 2)), tmp_vec, phinat_range_dofs);
    mv(M_pfield_, source_mu, tmp_vec, phinat_range_dofs);
    vector_axpy(range_phinat, 1. / Ca_, tmp_vec, phinat_range_dofs);
    // apply missing parts of A
    mv(M_nonlin_pfield_, source_phi, tmp_vec, mu_range_dofs);
    scal(tmp_vec, 1. / epsilon_, mu_range_dofs);
    add(range_mu, tmp_vec, mu_range_dofs);
    // apply G
    mv(Dphi_f_pfield_, source_phi, tmp_vec, phinat_range_dofs);
    add(range_phinat, tmp_vec, phinat_range_dofs);
  } else {
    fill_tmp_pfield(cell, source, restricted);
    assemble_nonlinear_part_of_pfield_residual(residual, cell, restricted);
  }
  if (restricted) {
    const auto& range_dofs = *pfield_deim_range_dofs_[cell];
    VectorType ret(range_dofs.size());
    for (size_t ii = 0; ii < range_dofs.size(); ++ii)
      ret[ii] = residual[range_dofs[ii]];
    return ret;
  } else {
    return residual;
  }
}

void CellModelSolver::update_pfield_parameters(
    const double Be, const double Ca, const double Pa, const size_t cell, const bool restricted)
{
  DUNE_THROW_IF(
      cell != 0, Dune::NotImplemented, "This may not work for several cells, check before removing this line!");
  if (XT::Common::FloatCmp::ne(Be, Be_) || XT::Common::FloatCmp::ne(Ca, Ca_)
      || XT::Common::FloatCmp::ne(Pa, last_pfield_Pa_)) {
    std::cout << "Pfield params updated, old (Be, Ca, Pa) = " << Be_ << ", " << Ca_ << ", " << last_pfield_Pa_
              << ", new = " << Be << ", " << Ca << ", " << Pa << std::endl;
    Be_ = Be;
    Ca_ = Ca;
    Pa_ = Pa;
    last_pfield_Pa_ = Pa_;
    XT::Common::Parameter param({{"gamma", {gamma_}}, {"epsilon", {epsilon_}}, {"Be", {Be_}}, {"Ca", {Ca_}}});
    pfield_jac_linear_op_.set_params(param, restricted);
    pfield_solver_.set_params(param, restricted);
  }
}

//******************************************************************************************************************
//******************************************* Apply inverse operators **********************************************
//******************************************************************************************************************

// Applies inverse stokes operator (solves F(y) = 0)
CellModelSolver::VectorType CellModelSolver::apply_inverse_stokes_operator() const
{
  EigenVectorType ret(size_u_ + size_p_);
  ret.backend() = stokes_solver_->solve(stokes_rhs_vector_.backend());
  return XT::Common::convert_to<VectorType>(ret);
}

// Applies inverse orientation field operator (solves F(y) = 0)
// y_guess is the initial guess for the Newton iteration
CellModelSolver::VectorType CellModelSolver::apply_inverse_ofield_operator(const VectorType& y_guess, const size_t cell)
{
  // *********** Newton ******************************
  const auto tol = 1e-12;
  const auto max_iter = 200;
  const auto max_dampening_iter = 1000;

  // auto l2_norm_P = l2_norm(grid_view_, P_[cell]);
  // auto l2_norm_Pnat = l2_norm(grid_view_, Pnat_[cell]);

  // ********* compute residual *********
  auto begin = std::chrono::steady_clock::now();
  auto residual = apply_ofield_operator(y_guess, cell);
  auto res_norm = ofield_residual_norm(residual);
  std::cout << "Newton: Initial Residual: " << res_norm << std::endl;
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

    // backtracking line search
    const double gamma = 0.001;
    while (candidate_res > (1 - gamma * lambda) * res_norm) {
      DUNE_THROW_IF(k >= max_dampening_iter,
                    Exceptions::operator_error,
                    "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                        << res_norm << "\nl = " << iter << "\n");
      y_n_plus_1 = y_n + update * lambda;
      residual = apply_ofield_operator(y_n_plus_1, cell);
      candidate_res = ofield_residual_norm(residual);
      lambda /= 2;
      k += 1;
    }
    y_n = y_n_plus_1;
    res_norm = candidate_res;
    std::cout << "Newton: Iter " << iter << " Current residual: " << res_norm << std::endl;
    iter += 1;
  } // while (true)
  return y_n;
}

// Applies inverse phase field operator (solves F(y) = 0)
// y_guess is the initial guess for the Newton iteration
CellModelSolver::VectorType CellModelSolver::apply_inverse_pfield_operator(const VectorType& y_guess, const size_t cell)
{
  // *********** Newton ******************************
  const auto tol = 1e-12;
  const auto max_iter = 200;
  const auto max_dampening_iter = 1000;

  // const auto l2_norm_phi = l2_norm(grid_view_, phi_[cell]);
  // const auto l2_norm_phinat = l2_norm(grid_view_, phinat_[cell]);
  // const auto l2_norm_mu = l2_norm(grid_view_, mu_[cell]);

  // ********* compute residual *********
  auto begin = std::chrono::steady_clock::now();
  auto residual = apply_pfield_operator(y_guess, cell, false);
  auto res_norm = pfield_residual_norm(residual);
  std::cout << "Newton: Initial Residual: " << res_norm << std::endl;
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

    if (iter > 10) {
      // apply damping via backtracking line search
      size_t k = 0;
      auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
      double lambda = 1;
      const double gamma = 0.001;
      while (candidate_res > (1 - gamma * lambda) * res_norm) {
        DUNE_THROW_IF(k >= max_dampening_iter,
                      Exceptions::operator_error,
                      "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                          << res_norm << "\nl = " << iter << "\n");
        x_n_plus_1 = x_n + update * lambda;
        residual = apply_pfield_operator(x_n_plus_1, cell, false);
        candidate_res = pfield_residual_norm(residual);
        // std::cout << "Candidate res: " << candidate_res << std::endl;
        lambda /= 2;
        k += 1;
      }
      x_n = x_n_plus_1;
      res_norm = candidate_res;
    } else {
      x_n += update;
      residual = apply_pfield_operator(x_n, cell, false);
      res_norm = pfield_residual_norm(residual);
    }
    std::cout << "Newton: Iter " << iter << " Current residual: " << res_norm << std::endl;
    iter += 1;
  } // while (true)
  return x_n;
}

//******************************************************************************************************************
//********************************************** Apply jacobians ***************************************************
//******************************************************************************************************************

void CellModelSolver::set_pfield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted)
{
  assemble_pfield_nonlinear_jacobian(source, cell, restricted);
}

void CellModelSolver::set_pfield_jacobian_state_dofs(const std::vector<R>& source, const size_t cell)
{
  const size_t expected_size = pfield_deim_source_dofs_[cell][0].size();
  DUNE_THROW_IF(source.size() != expected_size,
                XT::Common::Exceptions::wrong_input_given,
                "Source has incorrect size " + XT::Common::to_string(source.size()) + ", should be "
                    + XT::Common::to_string(expected_size) + "!");
  assemble_pfield_nonlinear_jacobian(XT::Common::convert_to<VectorType>(source), cell, true);
}

CellModelSolver::VectorType
CellModelSolver::apply_pfield_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  return apply_pfield_helper(y, cell, restricted, true);
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_pfield_jacobian(const VectorType& rhs, const size_t cell)
{
  return pfield_solver_.apply(rhs, cell);
}

void CellModelSolver::set_ofield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted)
{
  assemble_ofield_nonlinear_jacobian(source, cell, restricted);
}

void CellModelSolver::set_ofield_jacobian_state_dofs(const std::vector<R>& source, const size_t cell)
{
  const size_t expected_size = ofield_deim_source_dofs_[cell][1].size();
  DUNE_THROW_IF(source.size() != expected_size,
                XT::Common::Exceptions::wrong_input_given,
                "Source has incorrect size " + XT::Common::to_string(source.size()) + ", should be "
                    + XT::Common::to_string(expected_size) + "!");
  assemble_ofield_nonlinear_jacobian(XT::Common::convert_to<VectorType>(source), cell, true);
}

CellModelSolver::VectorType
CellModelSolver::apply_ofield_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  return apply_ofield_helper(y, cell, restricted, true);
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_ofield_jacobian(const VectorType& rhs, const size_t cell)
{
  return ofield_solver_.apply(rhs, cell);
}

CellModelSolver::VectorType CellModelSolver::apply_stokes_jacobian(const VectorType& y, const bool restricted)
{
  return apply_stokes_helper(y, restricted, true);
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
void CellModelSolver::assemble_stokes_rhs(const bool restricted)
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
        using TestBasisType = typename LocalElementIntegralFunctional<E, d>::GenericIntegrand::LocalTestBasisType;
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
          const auto Fa_inv_phi_tilde = -this->Fa_inv_ * phi_tilde;
          const auto xi_p_1 = 0.5 * (this->xi_ + 1);
          const auto xi_m_1 = 0.5 * (this->xi_ - 1);
          for (size_t ii = 0; ii < sz; ++ii) {
            for (size_t mm = 0; mm < d; ++mm) {
              const auto factor1_mm = Fa_inv_phi_tilde * P[mm] - xi_p_1 * Pnat[mm];
              const auto factor2_mm = xi_m_1 * P[mm];
              result[ii] += (phinat_grad_phi[mm] + grad_P_T_times_Pnat[mm]) * test_basis_values_[ii][mm];
              result[ii] += factor1_mm * (P * test_basis_grads_[ii][mm]);
              result[ii] -= factor2_mm * (Pnat * test_basis_grads_[ii][mm]);
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
  const auto scal = scal_func<EigenVectorViewType>(restricted);
  const auto& u_range_dofs = u_deim_range_dofs_;
  scal(stokes_f_vector_, 0., u_range_dofs);
  if (!restricted)
    f_functional.assemble(use_tbb_);
  else
    f_functional.assemble_range(stokes_deim_entities_);
  // apply dirichlet constraints for u.
  // TODO: add restricted version
  u_dirichlet_constraints_.apply(stokes_f_vector_);
}

// Computes orientation field rhs using currently stored values of variables and stores in ofield_rhs_vector_
void CellModelSolver::assemble_ofield_rhs(const size_t cell, const bool restricted)
{
  const auto& P_range_dofs = P_deim_range_dofs_[cell];
  const auto mv = mv_func<VectorViewType, VectorViewType>(restricted);
  mv(M_ofield_, P_[cell].dofs().vector(), ofield_f_vector_, P_range_dofs);
  if (XT::Common::FloatCmp::ne(beta_, 0.)) {
    const auto& Pnat_range_dofs = Pnat_deim_range_dofs_[cell];
    auto g_functional = make_vector_functional(P_space_, ofield_g_vector_);
    const auto scal = scal_func<VectorViewType>(restricted);
    scal(ofield_g_vector_, 0., Pnat_range_dofs);
    XT::Functions::GenericGridFunction<E, d> g(
        /*order = */ phi_space_.max_polorder(),
        /*post_bind_func*/
        [cell, this](const E& element) { this->bind_phi(cell, element); },
        /*evaluate_func*/
        [cell, factor = beta_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate rhs terms
          const auto grad_phi = this->grad_phi(cell, x_local, param);
          auto ret = grad_phi;
          ret *= factor;
          return ret;
        });
    g_functional.append(
        LocalElementIntegralFunctional<E, d>(LocalElementProductScalarWeightIntegrand<E, d>().with_ansatz(g)));
    if (!restricted)
      g_functional.assemble(use_tbb_);
    else
      g_functional.assemble_range(ofield_deim_entities_[cell]);
  }
}

// Computes phase field rhs using currently stored values of variables and stores in pfield_rhs_vector_
void CellModelSolver::assemble_pfield_rhs(const size_t cell, const bool restricted)
{
  const auto& phi_range_dofs = phi_deim_range_dofs_[cell];
  const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
  auto r1_functional = make_vector_functional(phi_space_, pfield_r1_vector_);
  // calculate r0
  const auto mv = mv_func<VectorType, VectorViewType>(restricted);
  auto& phi_n = phi_tmp_vec2_;
  phi_n = phi_[cell].dofs().vector();
  mv(M_pfield_, phi_n, pfield_r0_vector_, phi_range_dofs);

  // calculate h
  const auto scal = scal_func<VectorViewType>(restricted);
  scal(pfield_r1_vector_, 0., phinat_range_dofs);
  XT::Functions::GenericGridFunction<E, 1, 1> r1_pf(
      /*order = */ 2 * P_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_P(cell, element); },
      /*evaluate_func*/
      [cell, factor1 = -c_1_ / (2. * Pa_), factor2 = -beta_ / Pa_, this](const DomainType& x_local,
                                                                         const XT::Common::Parameter& param) {
        // evaluate P, divP
        const auto Pn = this->eval_P(cell, x_local, param);
        const auto grad_P = this->grad_P(cell, x_local, param);
        R div_P(0.);
        for (size_t ii = 0; ii < d; ++ii)
          div_P += grad_P[ii][ii];
        auto ret = factor1 * (Pn * Pn) + factor2 * div_P;
        return ret;
      });
  r1_functional.append(
      LocalElementIntegralFunctional<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(1.).with_ansatz(r1_pf)));
  // assemble rhs
  if (!restricted)
    r1_functional.assemble(use_tbb_);
  else
    r1_functional.assemble_range(pfield_deim_entities_[cell]);
}

// assembles linear part of orientation field jacobian and stores in S_ofield_
void CellModelSolver::assemble_ofield_linear_jacobian(const size_t cell, const bool restricted)
{
  const auto& P_range_dofs = P_deim_range_dofs_[cell];
  // calculate A
  // Omega - xi D = (1-xi)/2 \nabla u^T - (1+xi)/2 \nabla u
  XT::Functions::GenericGridFunction<E, d, d> Omega_minus_xi_D_transposed(
      /*order = */ u_space_.max_polorder(),
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
  set_mat_to_zero(A_ofield_, restricted, ofield_submatrix_pattern_, P_range_dofs);
  A_ofield_op_->clear();
  A_ofield_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(Omega_minus_xi_D_transposed)));
  A_ofield_op_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(u_)));
  C_ofield_linear_part_op_->clear();
  C_ofield_linear_part_op_->append(*A_ofield_op_);

  // calculate linear part S_10 = C
  // TODO: add restricted version
  C_ofield_linear_part_.backend() = C_ofield_elliptic_part_.backend();
  XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_inv_phi(
      /*order = */ phi_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_phi(cell, element); },
      /*evaluate_func*/
      [cell, factor = c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi = this->eval_phi(cell, x_local, param);
        return factor * phi;
      });
  C_ofield_linear_part_op_->append(
      LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(c1_Pa_inv_phi)));
  if (!restricted)
    C_ofield_linear_part_op_->assemble(use_tbb_);
  else
    C_ofield_linear_part_op_->assemble_range(ofield_deim_entities_[cell]);
  // In the restricted case, we cannot use the solver anyway
  if (ofield_solver_.is_schur_solver() && !restricted) {
    S_schur_ofield_linear_part_.backend() = M_ofield_.backend();
    S_schur_ofield_linear_part_.axpy(dt_, A_ofield_);
    S_schur_ofield_linear_part_.axpy(-dt_ / kappa_, C_ofield_linear_part_);
  }
}

// assembles nonlinear part of orientation field jacobian and adds to S_ofield_
// if assemble_ofield_linear_jacobian has been called first, S_ofield now contains the whole orientation field
// jacobian
void CellModelSolver::assemble_ofield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  fill_tmp_ofield(cell, y, restricted);
  assemble_C_ofield_nonlinear_part(cell, restricted);
  if (!restricted) {
    if (ofield_solver_.is_schur_solver()) {
      ofield_solver_.schur_matrix() = S_schur_ofield_linear_part_;
      ofield_solver_.schur_matrix().axpy(-dt_ / kappa_, C_ofield_nonlinear_part_);
    }
    ofield_solver_.prepare(cell, restricted);
  }
}

void CellModelSolver::assemble_C_ofield_nonlinear_part(const size_t cell, const bool restricted)
{
  XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_P2(
      /*order = */ 2. * P_space_.max_polorder(),
      /*post_bind_func*/
      [cell, this](const E& element) { this->bind_P(cell, element); },
      /*evaluate_func*/
      [cell, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto P_n = this->eval_P(cell, x_local, param);
        return factor * P_n.two_norm2();
      });
  set_mat_to_zero(C_ofield_nonlinear_part_, restricted, ofield_submatrix_pattern_, Pnat_deim_range_dofs_[cell]);
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
void CellModelSolver::assemble_pfield_linear_jacobian(const size_t cell, const bool restricted)
{
  // assemble matrix S_{00} = M + dt D
  assemble_B_pfield(cell, restricted);
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

void CellModelSolver::assemble_B_pfield(const size_t cell, const bool restricted)
{
  const auto& phi_range_dofs = phi_deim_range_dofs_[cell];
  set_mat_to_zero(B_pfield_, restricted, pfield_submatrix_pattern_, phi_range_dofs);
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
  B_pfield_op_->clear();
  B_pfield_op_->append(
      LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u)));
  if (!restricted)
    B_pfield_op_->assemble(use_tbb_);
  else
    B_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
}

// assembles nonlinear part of phase field jacobian
void CellModelSolver::assemble_pfield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  fill_tmp_pfield(cell, y, restricted);
  assemble_M_nonlin_pfield(cell, restricted);
  assemble_Dphi_f_pfield(cell, restricted);
}

// stores matrix with entries \int (3 phi^2 - 1) varphi_i varphi_j in M_nonlin_pfield_
void CellModelSolver::assemble_M_nonlin_pfield(const size_t cell, const bool restricted)
{
  const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
  const auto& mu_range_dofs = mu_deim_range_dofs_[cell];
  set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, mu_range_dofs);
  if (restricted)
    set_mat_to_zero(M_nonlin_pfield_, restricted, pfield_submatrix_pattern_, phinat_range_dofs);
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

// stores nonlinear part of block G of the phase field jacobian matrix in Dphi_f_pfield_nonlinear_part_
void CellModelSolver::assemble_Dphi_f_pfield(const size_t cell, const bool restricted)
{
  const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
  set_mat_to_zero(Dphi_f_pfield_, restricted, pfield_submatrix_pattern_, phinat_range_dofs);
  XT::Functions::GenericGridFunction<E, 1, 1> Dphi_f_prefactor(
      /*order = */ num_cells_ > 1 ? 2 * phi_space_.max_polorder() + 20 : 2 * phi_space_.max_polorder(),
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
  Dphi_f_pfield_op_->clear();
  Dphi_f_pfield_op_->append(
      LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(Dphi_f_prefactor)));
  if (!restricted)
    Dphi_f_pfield_op_->assemble(use_tbb_);
  else
    Dphi_f_pfield_op_->assemble_range(pfield_deim_entities_[cell]);
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
      /*order = */ num_cells_ > 1 ? 3 * phi_space_.max_polorder() + 20 : 3 * phi_space_.max_polorder(),
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
      LocalElementProductScalarWeightIntegrand<E, 1>().with_ansatz(nonlinear_res_pf1)));
  nonlinear_res2_functional.append(LocalElementIntegralFunctional<E, 1>(
      LocalElementProductScalarWeightIntegrand<E, 1>().with_ansatz(nonlinear_res_pf2)));
  nonlinear_res1_functional.append(nonlinear_res2_functional);
  if (!restricted)
    nonlinear_res1_functional.assemble(use_tbb_);
  else
    nonlinear_res1_functional.walk_range(pfield_deim_entities_[cell]);
}

//******************************************************************************************************************
//******************************************* DEIM related methods *************************************************
//******************************************************************************************************************

const std::vector<std::vector<size_t>>& CellModelSolver::stokes_deim_source_dofs() const
{
  return stokes_deim_source_dofs_;
}

const std::vector<std::vector<size_t>>& CellModelSolver::ofield_deim_source_dofs(const size_t cell) const
{
  return ofield_deim_source_dofs_[cell];
}

const std::vector<std::vector<size_t>>& CellModelSolver::pfield_deim_source_dofs(const size_t cell) const
{
  return pfield_deim_source_dofs_[cell];
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
  if (testcase == "single_cell" || testcase == "two_cells" || testcase == "channel")
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
    return {{30., 30.}};
  else if (testcase == "channel")
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
  if (testcase == "single_cell" || testcase == "channel")
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
  if (testcase == "single_cell" || testcase == "channel")
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

// sets temporary orientation field discrete functions to source values
void CellModelSolver::fill_tmp_ofield(const size_t cell, const VectorType& source, const bool restricted) const
{
  if (!restricted) {
    ConstVectorViewType P_vec(source, 0, size_P_);
    P_tmp_[cell].dofs().vector() = P_vec;
  } else {
    const auto& source_dofs = ofield_deim_source_dofs_[cell][1];
    const size_t Pnat_begin = Pnat_deim_source_dofs_begin_[cell];
    auto& P_vec = P_tmp_[cell].dofs().vector();
    if (source.size() == 2 * size_P_) {
      for (size_t ii = 0; ii < Pnat_begin; ++ii)
        P_vec.set_entry(source_dofs[ii], source[source_dofs[ii]]);
    } else {
      DUNE_THROW_IF(
          source.size() != source_dofs.size(), XT::Common::Exceptions::wrong_input_given, "Source has wrong size!");
      for (size_t ii = 0; ii < Pnat_begin; ++ii)
        P_vec.set_entry(source_dofs[ii], source[ii]);
    }
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
    const auto& source_dofs = pfield_deim_source_dofs_[cell][0];
    const size_t phinat_begin = phinat_deim_source_dofs_begin_[cell];
    const size_t mu_begin = mu_deim_source_dofs_begin_[cell];
    auto& phi_vec = phi_tmp_[cell].dofs().vector();
    auto& mu_vec = mu_tmp_[cell].dofs().vector();
    if (source.size() == 3 * size_phi_) {
      for (size_t ii = 0; ii < phinat_begin; ++ii)
        phi_vec.set_entry(source_dofs[ii], source[source_dofs[ii]]);
      for (size_t ii = mu_begin; ii < source_dofs.size(); ++ii)
        mu_vec.set_entry(source_dofs[ii] - 2 * size_phi_, source[source_dofs[ii]]);
    } else {
      DUNE_THROW_IF(
          source.size() != source_dofs.size(), XT::Common::Exceptions::wrong_input_given, "Source has wrong size!");
      for (size_t ii = 0; ii < phinat_deim_source_dofs_begin_[cell]; ++ii)
        phi_vec.set_entry(source_dofs[ii], source[ii]);
      for (size_t ii = mu_deim_source_dofs_begin_[cell]; ii < source_dofs.size(); ++ii)
        mu_vec.set_entry(source_dofs[ii] - 2 * size_phi_, source[ii]);
    }
  }
}

// error norm used in orientation field Newton iteration
// TODO: use appropriate norm
double CellModelSolver::ofield_residual_norm(const VectorType& residual) const
{
  return residual.l2_norm();
}

// error norm used in phase field Newton iteration
// TODO: use appropriate norm
double CellModelSolver::pfield_residual_norm(const VectorType& residual) const
{
  return residual.l2_norm();
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
