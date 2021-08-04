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
#include <dune/xt/la/solver/istl/saddlepoint.hh>

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

#include <boost/geometry.hpp>

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
                                 const double gamma, // phase field mobility coefficient
                                 const double epsilon, // phase field parameter
                                 const int overintegrate,
                                 const CellModelLinearSolverType pfield_solver_type,
                                 const CellModelMassMatrixSolverType pfield_mass_matrix_solver_type,
                                 const CellModelLinearSolverType ofield_solver_type,
                                 const CellModelMassMatrixSolverType ofield_mass_matrix_solver_type,
                                 const StokesSolverType stokes_solver_type,
                                 const double outer_reduction,
                                 const int outer_restart,
                                 const int outer_verbose,
                                 const double inner_reduction,
                                 const int inner_maxit,
                                 const int inner_verbose,
                                 const double bending,
                                 const double conservative)
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
  , gamma_(gamma)
  , Be_(Be)
  , Ca_(Ca)
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
  , num_mutexes_u_(use_tbb ? size_u_ / 100 : 0)
  , num_mutexes_ofield_(use_tbb ? size_P_ / 100 : 0)
  , num_mutexes_pfield_(use_tbb ? size_phi_ / 100 : 0)
  , num_pfield_variables_(!bending && !conservative ? 1 : 3) // TODO: other combinations of bending and conservative
  , stokes_vector_(size_u_ + size_p_, 0.)
  , ofield_vectors_(num_cells_, VectorType(2 * size_P_, 0.))
  , pfield_vectors_(num_cells_, VectorType(num_pfield_variables_ * size_phi_, 0.))
  , u_view_(stokes_vector_, 0, size_u_)
  , p_view_(stokes_vector_, size_u_, size_u_ + size_p_)
  , u_(u_space_, u_view_, "u")
  , p_(p_space_, p_view_, "p")
  , A_stokes_(size_u_, size_u_, make_element_sparsity_pattern(u_space_, u_space_, grid_view_), 100)
  , BT_stokes_(size_u_, size_p_, make_element_sparsity_pattern(u_space_, p_space_, grid_view_), 100)
  , stokes_solver_type_(stokes_solver_type)
  , stokes_solver_(std::make_shared<LUSolverType>())
  , stokes_A_solver_(std::make_shared<LDLTSolverType>())
  , stokes_jac_linear_op_(*this)
  , stokes_rhs_vector_(size_u_ + size_p_, 0., num_mutexes_u_)
  , stokes_f_vector_(stokes_rhs_vector_, 0, size_u_)
  , stokes_g_vector_(stokes_rhs_vector_, size_u_, size_u_ + size_p_)
  // , p_basis_integrated_vector_(size_p_)
  , stokes_tmp_vec_(size_u_ + size_p_)
  , stokes_tmp_vec2_(size_u_ + size_p_)
  , stokes_p_tmp_vec_(size_p_)
  , stokes_p_tmp_vec2_(size_p_)
  , stokes_u_tmp_vec_(size_u_)
  , stokes_u_tmp_vec2_(size_u_)
  , u_dirichlet_constraints_(make_dirichlet_constraints(u_space_, boundary_info_))
  , ofield_submatrix_pattern_(make_element_sparsity_pattern(P_space_, P_space_, grid_view_))
  , M_ofield_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , B_ofield_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , E_ofield_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , C_ofield_incl_coeffs_and_sign_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , Dd_f_ofield_incl_coeffs_and_sign_(size_P_, size_P_, ofield_submatrix_pattern_, num_mutexes_ofield_)
  , S_schur_ofield_linear_part_(size_P_, size_P_, ofield_submatrix_pattern_, 0)
  , ofield_jac_linear_op_(*this)
  , ofield_solver_(dt_,
                   kappa_,
                   Pa_,
                   M_ofield_,
                   E_ofield_,
                   B_ofield_,
                   C_ofield_incl_coeffs_and_sign_,
                   Dd_f_ofield_incl_coeffs_and_sign_,
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
  , ofield_tmp_vec_(2 * size_P_, 0.)
  , ofield_tmp_vec2_(2 * size_P_, 0., num_mutexes_ofield_)
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
  , E_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , Dmu_f_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , Dphi_f_pfield_incl_coeffs_and_sign_(size_phi_, size_phi_, pfield_submatrix_pattern_, num_mutexes_pfield_)
  , pfield_jac_linear_op_(*this)
  , pfield_solver_(num_pfield_variables_,
                   dt_,
                   gamma_,
                   epsilon_,
                   Be_,
                   Ca_,
                   M_pfield_,
                   E_pfield_,
                   B_pfield_,
                   Dphi_f_pfield_incl_coeffs_and_sign_,
                   Dmu_f_pfield_,
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
  , pfield_rhs_vector_(num_pfield_variables_ * size_phi_, 0., num_mutexes_pfield_)
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
  , pfield_tmp_vec_(num_pfield_variables_ * size_phi_, 0.)
  , pfield_tmp_vec2_(num_pfield_variables_ * size_phi_, 0., num_mutexes_pfield_)
  , phi_tmp_vec_(size_phi_, 0.)
  , phi_tmp_vec2_(size_phi_, 0.)
  , u_tmp_vec_(size_u_, 0.)
  , P_tmp_vec_(size_P_, 0.)
  , u_tmp_(u_space_)
  , partitioning_(
        grid_view_,
        use_tbb ? DXTC_CONFIG_GET("threading.partition_factor", 1u) * XT::Common::threadManager().current_threads() : 1)
  , bending_(bending)
  , conservative_(conservative)
{
  std::cout << dx_ << ", " << epsilon_ << std::endl;
  std::cout << "DOFs: Pfield " << XT::Common::to_string(3 * size_phi_) << ", Ofield "
            << XT::Common::to_string(2 * size_P_) << ", Stokes " << XT::Common::to_string(size_u_ + size_p_)
            << std::endl;

  /************************** create and project initial values*****************************************
   ************************** we only need initial values for P and phi *******************************/

  std::shared_ptr<const XT::Functions::FunctionInterface<d, d>> u_initial_func;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d>>> phi_initial_funcs;
  std::vector<std::shared_ptr<const XT::Functions::FunctionInterface<d, d>>> P_initial_funcs;

  // interpolate initial and boundary values
  if (testcase == "single_cell" || testcase == "channel" || testcase == "single_cell_dirichlet") {
    // Initially, cell is circular with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 30]^2
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(15, 15) - x|.
    FieldVector<double, d> center{upper_right_[0] / 2., upper_right_[1] / 2.};
    // FieldVector<double, d> center{0., upper_right_[1] / 2.};
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
    //     P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
    //         50,
    //         /*evaluate=*/
    //         [&phi_in = phi_initial_funcs[0]](const auto& x, const auto& param) {
    //           // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
    //           // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
    //           // auto ret = FieldVector<double, d>({1. + rand1, 0. +
    //           // rand2});
    //           auto ret = FieldVector<double, d>({1., 0.});
    //           ret *= (phi_in->evaluate(x, param) + 1.) / 2.;
    //           return ret;
    //         },
    //         /*name=*/"P_initial"));
    P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
        50,
        /*evaluate=*/
        [&phi_in = phi_initial_funcs[0]](const auto& x, const auto& param) {
          // auto rand1 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto rand2 = ((std::rand() % 2000) - 1000) / 20000.;
          // auto ret = FieldVector<double, d>({1. + rand1, 0. +
          // rand2});
          // double rand1 = std::rand() - RAND_MAX / 2;
          // double rand2 = std::rand() - RAND_MAX / 2;
          // rand1 /= std::sqrt(rand1 * rand1 + rand2 * rand2);
          // rand2 /= std::sqrt(rand1 * rand1 + rand2 * rand2);
          // auto ret = FieldVector<double, d>({rand1, rand2});
          // small bias to point more in x direction
          // if (((x[0] > 9.5 && x[0] < 10.5) || (x[0] > 19.5 && x[0] < 20.5)) && x[1] > 14.5 && x[1] < 15.5)
          // ret[0] += (ret[0] > 0) ? 0.05 : -0.05;
          auto ret = FieldVector<double, d>({1., 0.});
          ret *= (phi_in->evaluate(x, param) + 1.) / 2.;
          return ret;
        },
        /*name=*/"P_initial"));
    u_initial_func = std::make_shared<const XT::Functions::ConstantFunction<d, d>>(0.);
  } else if (testcase == "cell_isolation_experiment") {
    // Initially, cell is elongated with Radius R=5 and placed in the center of the domain
    // \Omega = [0, 30]^2
    // Initial condition for \phi thus is \tanh(\frac{r}{\sqrt{2}\epsilon}) with r the signed distance function to the
    // membrane, i.e. r(x) = 5 - |(15, 15) - x|.
    FieldVector<double, d> center{upper_right_[0] / 2., upper_right_[1] / 2.};
    using BoostPointType = boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian>;
    using BoostRingType = boost::geometry::model::ring<BoostPointType>;
    using BoostPolygonType = boost::geometry::model::polygon<BoostPointType>;
    // A ring does not have holes
    BoostRingType ring{{15., 25.}, {15., 20.}, {25., 15.}, {27., 20.}, {26., 22.}, {15., 25.}};
    // BoostRingType ring{{15.2, 25.}, {15., 24.8}, {15., 20.}, {25., 15.}, {27., 20.}, {26., 22.}, {15.2, 25.}};
    // A polygon can have holes, here we add the same inner and outer ring to only get the boundary of the cell
    BoostPolygonType polygon{ring, ring};
    auto r = [ring, polygon](const auto& xr) {
      BoostPointType point{xr[0], xr[1]};
      const double distance = boost::geometry::distance(polygon, point);
      return boost::geometry::within(point, ring) ? distance : -distance;
    };
    phi_initial_funcs.emplace_back(std::make_shared<XT::Functions::GenericFunction<d>>(
        50,
        /*evaluate=*/
        [r, epsilon = epsilon_](const auto& x, const auto& /*param*/) {
          return std::tanh(r(x) / (std::sqrt(2.) * epsilon));
        },
        /*name=*/"phi_initial"));

    P_initial_funcs.emplace_back(std::make_shared<const XT::Functions::GenericFunction<d, d>>(
        50,
        /*evaluate=*/
        [&phi_in = phi_initial_funcs[0]](const auto& x, const auto& param) {
          // auto ret = FieldVector<double, d>({1., 0.});
          auto ret = FieldVector<double, d>({-0.82, 0.57});
          ret *= (phi_in->evaluate(x, param) + 1.) / 2.;
          return ret;
        },
        /*name=*/"P_initial"));
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
    if (num_pfield_variables_ > 1) {
      phinat_view_.emplace_back(pfield_vectors_[kk], size_phi_, 2 * size_phi_);
      mu_view_.emplace_back(pfield_vectors_[kk], 2 * size_phi_, 3 * size_phi_);
    }
    const auto kk_str = XT::Common::to_string(kk);
    P_.emplace_back(make_discrete_function(P_space_, P_view_[kk], "P_" + kk_str));
    Pnat_.emplace_back(make_discrete_function(P_space_, Pnat_view_[kk], "Pnat_" + kk_str));
    phi_.emplace_back(make_discrete_function(phi_space_, phi_view_[kk], "phi_" + kk_str));
    if (num_pfield_variables_ > 1) {
      phinat_.emplace_back(make_discrete_function(phi_space_, phinat_view_[kk], "phinat_" + kk_str));
      mu_.emplace_back(make_discrete_function(phi_space_, mu_view_[kk], "mu_" + kk_str));
    }
    P_tmp_.emplace_back(P_space_);
    Pnat_tmp_.emplace_back(P_space_);
    phi_tmp_.emplace_back(phi_space_);
    if (num_pfield_variables_ > 1) {
      phinat_tmp_.emplace_back(phi_space_);
      mu_tmp_.emplace_back(phi_space_);
    }
    default_interpolation(*phi_initial_funcs[kk], phi_[kk]);
    default_interpolation(*P_initial_funcs[kk], P_[kk]);
  }

  /*************************************************************************************************
   *************************************** Stokes **************************************************
   *************************************************************************************************/


  MatrixOperator<MatrixType, PGV, d> A_stokes_op(grid_view_, u_space_, u_space_, A_stokes_);
  MatrixOperator<MatrixType, PGV, 1, 1, d> BT_stokes_op(grid_view_, p_space_, u_space_, BT_stokes_);
  MatrixType M_p_stokes(size_p_, size_p_, make_element_sparsity_pattern(p_space_, p_space_, grid_view_), 100);
  if (stokes_solver_type_ == StokesSolverType::schur_cg_A_direct_prec_mass
      || stokes_solver_type_ == StokesSolverType::schur_cg_A_direct_prec_masslumped) {
    MatrixOperator<MatrixType, PGV, 1, 1, 1> M_p_stokes_op(grid_view_, p_space_, p_space_, M_p_stokes);
    M_p_stokes_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>()));
    A_stokes_op.append(M_p_stokes_op);
  }
  // Stokes matrix is [A B^T; B C]
  // calculate A_{ij} as \int 0.5 \nabla v_i \nabla v_j
  A_stokes_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(0.5)));
  // calculate (B^T)_{ij} as \int p_i div(v_j)
  BT_stokes_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
      LocalElementAnsatzValueTestDivProductIntegrand<E>()));

  // auto p_basis_integrated_functional = make_vector_functional(p_space_, p_basis_integrated_vector_);
  // const XT::Functions::ConstantGridFunction<E> one_function(1);
  // p_basis_integrated_functional.append(
  // LocalElementIntegralFunctional<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>().with_ansatz(one_function)));
  // BT_stokes_op.append(p_basis_integrated_functional);

  // Dirichlet constrainst for u
  A_stokes_op.append(u_dirichlet_constraints_);
  // assemble everything
  A_stokes_op.append(BT_stokes_op);
  std::cout << "Assembling stokes..." << std::flush;
  A_stokes_op.assemble(use_tbb_);
  std::cout << "done" << std::endl;


  // apply dirichlet constraints for u
  u_dirichlet_constraints_.apply(A_stokes_, /*clear_only = */ false, /* ensure_symmetry = */ true);
  for (const auto& DoF : u_dirichlet_constraints_.dirichlet_DoFs())
    BT_stokes_.clear_row(DoF);

  // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the u_size_-th row of
  // [A B^T; B C] to the unit vector. Matrix C is adapted below for the direct solver and when applying the Schur
  // complement for the schur solvers.
  stokes_g_vector_.set_entry(0, 0.);
  BT_stokes_.clear_col(0);
  if (stokes_solver_type_ == StokesSolverType::direct) {
    // Setup system matrix
    MatrixType S_stokes(size_u_ + size_p_, size_u_ + size_p_, create_stokes_pattern(u_space_, p_space_), 100);
    MatrixViewType A_stokes_view(S_stokes, 0, size_u_, 0, size_u_);
    MatrixViewType BT_stokes_view(S_stokes, 0, size_u_, size_u_, size_u_ + size_p_);
    MatrixViewType B_stokes_view(S_stokes, size_u_, size_u_ + size_p_, 0, size_u_);
    A_stokes_view = A_stokes_;
    BT_stokes_view = BT_stokes_;
    const auto BT_pattern = BT_stokes_.pattern();
    for (size_t ii = 0; ii < size_u_; ii++)
      for (const auto& jj : BT_pattern.inner(ii))
        B_stokes_view.set_entry(jj, ii, BT_stokes_.get_entry(ii, jj));
    S_stokes.set_entry(size_u_, size_u_, 1.); // C matrix
    // Factor system matrix
    std::cout << "Factorizing Stokes system matrix..." << std::flush;
    auto t1 = std::chrono::high_resolution_clock::now();
    stokes_solver_->compute(S_stokes.backend());
    stokes_solver_statistics_.setup_time_ = std::chrono::high_resolution_clock::now() - t1;
    if (stokes_solver_->info() != ::Eigen::Success)
      DUNE_THROW(Dune::InvalidStateException, "Failed to invert stokes system matrix!");
    std::cout << "done" << std::endl;
  } else {
    std::cout << "Factorizing stokes A matrix..." << std::flush;
    auto t1 = std::chrono::high_resolution_clock::now();
    stokes_A_solver_->compute(A_stokes_.backend());
    stokes_solver_statistics_.setup_time_ = std::chrono::high_resolution_clock::now() - t1;
    if (stokes_A_solver_->info() != ::Eigen::Success)
      DUNE_THROW(Dune::InvalidStateException, "Failed to invert stokes A matrix!");
    std::cout << "done" << std::endl;
    if (stokes_solver_type == StokesSolverType::schur_cg_A_direct) {
      stokes_preconditioner_ = std::make_shared<XT::LA::IdentityPreconditioner<StokesSchurComplementType>>(
          SolverCategory::Category::sequential);
    } else if (stokes_solver_type_ == StokesSolverType::schur_cg_A_direct_prec_mass) {
      std::cout << "Factorizing pressure mass matrix..." << std::flush;
      M_p_stokes_solver_ = std::make_shared<LDLTSolverType>();
      stokes_preconditioner_ = std::make_shared<StokesMassMatrixPreconditionerType>(M_p_stokes_solver_);
      t1 = std::chrono::high_resolution_clock::now();
      M_p_stokes_solver_->compute(M_p_stokes.backend());
      stokes_solver_statistics_.setup_time_ += std::chrono::high_resolution_clock::now() - t1;
      if (M_p_stokes_solver_->info() != ::Eigen::Success)
        DUNE_THROW(Dune::InvalidStateException, "Failed to invert M_p_stokes matrix!");
      std::cout << "done" << std::endl;
    } else if (stokes_solver_type_ == StokesSolverType::schur_cg_A_direct_prec_masslumped) {
      EigenVectorType M_p_stokes_lumped(size_p_, 0.);
      const auto M_p_pattern = make_element_sparsity_pattern(p_space_, p_space_, grid_view_);
      for (size_t ii = 0; ii < size_p_; ++ii)
        for (const size_t jj : M_p_pattern.inner(ii))
          M_p_stokes_lumped[ii] += M_p_stokes.get_entry(ii, jj);
      stokes_preconditioner_ = std::make_shared<LumpedMassMatrixPreconditioner<EigenVectorType>>(M_p_stokes_lumped);
    }
    stokes_schur_complement_ = std::make_shared<StokesSchurComplementType>(*this);
    stokes_schur_solver_ = std::make_shared<Dune::CGSolver<EigenVectorType>>(
        *stokes_schur_complement_, *stokes_preconditioner_, 1e-13, 1000, outer_verbose, false);
  }

  /*************************************************************************************************
   ************************************ Orientationfield *******************************************
   *************************************************************************************************/
  // calculate M_{ij} as \int \psi_i phi_j
  MatrixOperator<MatrixType, PGV, d> M_ofield_op(grid_view_, P_space_, P_space_, M_ofield_);
  M_ofield_op.append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductScalarWeightIntegrand<E, d>(1.)));
  E_ofield_ *= 0.;
  MatrixOperator<MatrixType, PGV, d> E_ofield_op(grid_view_, P_space_, P_space_, E_ofield_);
  E_ofield_op.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>()));
  M_ofield_op.append(E_ofield_op);

  std::cout << "Assembling ofield..." << std::flush;
  M_ofield_op.assemble(use_tbb_);
  std::cout << "done" << std::endl;
  std::cout << "Setting up ofield_solver..." << std::flush;
  auto t1 = std::chrono::high_resolution_clock::now();
  ofield_solver_.setup();
  ofield_solver_statistics_.setup_time_ = std::chrono::high_resolution_clock::now() - t1;
  std::cout << "done" << std::endl;

  /*************************************************************************************************
   **************************************** Phasefield *********************************************
   *************************************************************************************************/

  MatrixOperator<MatrixType, PGV, 1> M_pfield_op(grid_view_, phi_space_, phi_space_, M_pfield_);
  M_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductScalarWeightIntegrand<E, 1>(1.)));
  MatrixOperator<MatrixType, PGV, 1> E_pfield_op(grid_view_, phi_space_, phi_space_, E_pfield_);
  E_pfield_op.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
  M_pfield_op.append(E_pfield_op);
  std::cout << "Assembling pfield..." << std::flush;
  M_pfield_op.assemble(use_tbb_);
  std::cout << "done" << std::endl;
  std::cout << "Setting up pfield_solver..." << std::flush;
  t1 = std::chrono::high_resolution_clock::now();
  pfield_solver_.setup();
  pfield_solver_statistics_.setup_time_ = std::chrono::high_resolution_clock::now() - t1;
  std::cout << "done" << std::endl;


  /*************************************************************************************************
   ****************************************** Indices **********************************************
   *************************************************************************************************/

  global_indices_phi_ = std::vector<DynamicVector<size_t>>(grid_view_.size(0));
  global_indices_u_ = global_indices_phi_;
  global_indices_P_ = global_indices_phi_;
  for (const auto& element : Dune::elements(grid_view_)) {
    const auto index = grid_view_.indexSet().index(element);
    phi_space_.mapper().global_indices(element, global_indices_phi_[index]);
    u_space_.mapper().global_indices(element, global_indices_u_[index]);
    P_space_.mapper().global_indices(element, global_indices_P_[index]);
  }
  stokes_rhs_quad_ = std::make_shared<QuadratureStorage>(
      *this,
      std::max(phi_space_.max_polorder() + 2 * P_space_.max_polorder(), 2 * phi_space_.max_polorder())
          + u_space_.max_polorder() + overintegrate,
      true,
      true,
      true,
      true,
      false,
      true,
      true);
  pfield_rhs_quad_ =
      std::make_shared<QuadratureStorage>(*this,
                                          2 * P_space_.max_polorder() + phi_space_.max_polorder() + overintegrate,
                                          false,
                                          false,
                                          true,
                                          false,
                                          false,
                                          true,
                                          false);
  ofield_prepare_quad_ = std::make_shared<QuadratureStorage>(
      *this,
      std::max(u_space_.max_polorder(), phi_space_.max_polorder()) + 2 * P_space_.max_polorder() + overintegrate,
      true,
      true,
      true,
      true,
      false,
      true,
      false);
  ofield_residual_and_nonlinear_S10_quad_ = std::make_shared<QuadratureStorage>(
      *this, 4 * P_space_.max_polorder() + overintegrate, false, false, true, false, true, false, false);
  pfield_B_quad_ =
      std::make_shared<QuadratureStorage>(*this,
                                          u_space_.max_polorder() + 2 * phi_space_.max_polorder() + overintegrate,
                                          true,
                                          false,
                                          false,
                                          false,
                                          false,
                                          true,
                                          true);
  pfield_residual_and_nonlinear_jac_quad_ =
      std::make_shared<QuadratureStorage>(*this,
                                          4 * phi_space_.max_polorder() + overintegrate + 20 * (num_cells_ > 1),
                                          false,
                                          false,
                                          false,
                                          false,
                                          false,
                                          true,
                                          false);
} // constructor

// Computes the nonlinear part of the lower left block of the ofield matrix (i.e., -c_1/Pa D_d f(d))
struct CellModelSolver::OfieldNonlinearS10Functor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = OfieldNonlinearS10Functor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  OfieldNonlinearS10Functor(CellModelSolver& cellmodel, const QuadratureStorage& quad_storage)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , P_(0.)
    , P_times_P_basis_(quad_storage_.P_basis_size_, 0.)
    , result_(quad_storage_.P_basis_size_ * quad_storage_.P_basis_size_)
    , P_coeffs_(quad_storage_.P_basis_size_)
    , minus_frac_c1_Pa_(-cellmodel_.c_1_ / cellmodel_.Pa_)
  {}

  OfieldNonlinearS10Functor(const ThisType& other) = default;
  OfieldNonlinearS10Functor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new OfieldNonlinearS10Functor(*this);
  }

  void apply_local(const E& element) final
  {
    std::fill(result_.begin(), result_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_P = cellmodel_.global_indices_P_[index];
    const auto P_basis_size = quad_storage_.P_basis_size_;
    for (size_t ii = 0; ii < P_basis_size; ++ii)
      P_coeffs_[ii] = quad_storage_.P_dofs_[global_indices_P[ii]];
    for (size_t ll = 0; ll < quad_storage_.quadrature_.size(); ++ll) {
      const auto& quad_point = quad_storage_.quadrature_[ll];
      const auto factor = quad_storage_.integration_factors_[ll] * quad_point.weight();
      P_ *= 0.;
      for (size_t ii = 0; ii < P_basis_size; ++ii) {
        P_.axpy(P_coeffs_[ii], quad_storage_.P_basis_evaluated_[ll][ii]);
      }
      const auto factor1 = minus_frac_c1_Pa_ * P_.two_norm2() * factor;
      const auto factor2 = 2. * minus_frac_c1_Pa_ * factor;
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        P_times_P_basis_[ii] = P_ * quad_storage_.P_basis_evaluated_[ll][ii];
      for (size_t ii = 0; ii < P_basis_size; ++ii) {
        const auto offset_ii = ii * P_basis_size;
        const auto factor2_ii = P_times_P_basis_[ii] * factor2;
        for (size_t jj = 0; jj < P_basis_size; ++jj)
          result_[offset_ii + jj] +=
              quad_storage_.P_basis_products_[ll][ii][jj] * factor1 + P_times_P_basis_[jj] * factor2_ii;
      } // ii
    } // ll
    for (size_t ii = 0; ii < P_basis_size; ++ii) {
      const auto offset_ii = ii * P_basis_size;
      const auto index_ii = global_indices_P[ii];
      for (size_t jj = 0; jj < P_basis_size; ++jj)
        cellmodel_.Dd_f_ofield_incl_coeffs_and_sign_.add_to_entry(
            index_ii, global_indices_P[jj], result_[offset_ii + jj]);
    }
  }

  CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  DomainType P_;
  std::vector<double> P_times_P_basis_;
  std::vector<double> result_;
  std::vector<double> P_coeffs_;
  const double minus_frac_c1_Pa_;
}; // struct OfieldNonlinearS10Functor

// Computes the two ofield matrices B(u) and C(phi)
struct CellModelSolver::OfieldPrepareFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = OfieldPrepareFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  OfieldPrepareFunctor(CellModelSolver& cellmodel, const QuadratureStorage& quad_storage)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , u_(0.)
    , P_(0.)
    , grad_u_(0.)
    , J_inv_T_(0.)
    , u_basis_jacobians_transformed_(quad_storage_.u_basis_jacobians_)
    , P_basis_jacobians_transformed_(quad_storage_.P_basis_jacobians_)
    , P_basis_grad_times_u_(quad_storage_.P_basis_size_)
    , omega_minus_xi_D_T_basis_ii_(0.)
    , result1_(quad_storage_.P_basis_size_ * quad_storage_.P_basis_size_)
    , result2_(quad_storage_.P_basis_size_ * quad_storage_.P_basis_size_)
    , u_coeffs_(quad_storage_.u_basis_size_)
    , P_coeffs_(quad_storage_.P_basis_size_)
    , phi_coeffs_(quad_storage_.phi_basis_size_)
    , frac_c1_Pa_(cellmodel_.c_1_ / cellmodel_.Pa_)
    , grad_u_factor_((1. - cellmodel_.xi_) / 2.)
    , grad_u_T_factor_(-(1. + cellmodel_.xi_) / 2.)
  {}

  OfieldPrepareFunctor(const ThisType& other) = default;
  OfieldPrepareFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new OfieldPrepareFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_P = cellmodel_.global_indices_P_[index];
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const auto& global_indices_u = cellmodel_.global_indices_u_[index];
    const size_t u_basis_size = quad_storage_.u_basis_size_;
    const size_t P_basis_size = quad_storage_.P_basis_size_;
    const size_t phi_basis_size = quad_storage_.phi_basis_size_;
    const auto& quadrature = quad_storage_.quadrature_;
    for (size_t ii = 0; ii < u_basis_size; ++ii)
      u_coeffs_[ii] = quad_storage_.u_dofs_[global_indices_u[ii]];
    for (size_t ii = 0; ii < P_basis_size; ++ii)
      P_coeffs_[ii] = quad_storage_.P_dofs_[global_indices_P[ii]];
    for (size_t ii = 0; ii < phi_basis_size; ++ii)
      phi_coeffs_[ii] = quad_storage_.phi_dofs_[global_indices_phi[ii]];
    std::fill(result1_.begin(), result1_.end(), 0.);
    std::fill(result2_.begin(), result2_.end(), 0.);
    const auto& geometry = element.geometry();
    const bool geometry_affine = geometry.affine();
    J_inv_T_ = geometry.jacobianInverseTransposed(geometry.center());
    if (quad_storage_.P_affine_ && geometry_affine) {
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        for (size_t rr = 0; rr < d; ++rr)
          J_inv_T_.mv(quad_storage_.P_basis_jacobians_[0][ii][rr], P_basis_jacobians_transformed_[0][ii][rr]);
    }
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto factor = quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      if (!geometry_affine)
        J_inv_T_ = geometry.jacobianInverseTransposed(quadrature[ll].position());
      for (size_t ii = 0; ii < u_basis_size; ++ii)
        for (size_t rr = 0; rr < d; ++rr)
          J_inv_T_.mv(quad_storage_.u_basis_jacobians_[ll][ii][rr], u_basis_jacobians_transformed_[ll][ii][rr]);
      if (!quad_storage_.P_affine_ || !geometry_affine)
        for (size_t ii = 0; ii < P_basis_size; ++ii)
          for (size_t rr = 0; rr < d; ++rr)
            J_inv_T_.mv(quad_storage_.P_basis_jacobians_[ll][ii][rr], P_basis_jacobians_transformed_[ll][ii][rr]);

      u_ *= 0.;
      grad_u_ *= 0.;
      for (size_t ii = 0; ii < u_basis_size; ++ii) {
        u_.axpy(u_coeffs_[ii], quad_storage_.u_basis_evaluated_[ll][ii]);
        grad_u_.axpy(u_coeffs_[ii], u_basis_jacobians_transformed_[ll][ii]);
      }
      const double phi = XT::Common::transform_reduce(
          phi_coeffs_.begin(), phi_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
      const auto& P_basis_grads = quad_storage_.P_affine_ && geometry_affine ? P_basis_jacobians_transformed_[0]
                                                                             : P_basis_jacobians_transformed_[ll];
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        P_basis_grads[ii].mv(u_, P_basis_grad_times_u_[ii]);
      const auto factor2 = frac_c1_Pa_ * phi * factor;
      for (size_t ii = 0; ii < P_basis_size; ++ii) {
        const auto offset_ii = ii * P_basis_size;
        omega_minus_xi_D_T_basis_ii_ *= 0;
        // perform matrix vector multiplication
        for (size_t mm = 0; mm < d; ++mm)
          for (size_t nn = 0; nn < d; ++nn)
            omega_minus_xi_D_T_basis_ii_[mm] += (grad_u_[mm][nn] * grad_u_factor_ + grad_u_[nn][mm] * grad_u_T_factor_)
                                                * quad_storage_.P_basis_evaluated_[ll][ii][nn];
        tmp_vec1_ = omega_minus_xi_D_T_basis_ii_;
        tmp_vec1_ *= factor;
        tmp_vec2_ = quad_storage_.P_basis_evaluated_[ll][ii];
        tmp_vec2_ *= factor;
        for (size_t jj = 0; jj < P_basis_size; ++jj)
          result1_[offset_ii + jj] +=
              quad_storage_.P_basis_evaluated_[ll][jj] * tmp_vec1_ + P_basis_grad_times_u_[jj] * tmp_vec2_;
        const auto factor2_ii = quad_storage_.P_basis_evaluated_[ll][ii] * factor2;
        for (size_t jj = 0; jj < P_basis_size; ++jj)
          result2_[offset_ii + jj] += quad_storage_.P_basis_evaluated_[ll][jj] * factor2_ii;
      } // ii
    } // ll
    for (size_t ii = 0; ii < P_basis_size; ++ii) {
      const auto offset_ii = ii * P_basis_size;
      const auto index_ii = global_indices_P[ii];
      for (size_t jj = 0; jj < P_basis_size; ++jj) {
        const auto index_jj = global_indices_P[jj];
        cellmodel_.B_ofield_.add_to_entry(index_ii, index_jj, result1_[offset_ii + jj]);
        cellmodel_.C_ofield_incl_coeffs_and_sign_.add_to_entry(index_ii, index_jj, result2_[offset_ii + jj]);
      } // jj
    } // ii
  } // ... apply_local(...)

  CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  DomainType u_;
  DomainType P_;
  FieldMatrix<double, 2, 2> grad_u_;
  FieldMatrix<double, 2, 2> J_inv_T_;
  DomainType tmp_vec1_;
  DomainType tmp_vec2_;
  std::vector<std::vector<FieldMatrix<double, 2, 2>>> u_basis_jacobians_transformed_;
  std::vector<std::vector<FieldMatrix<double, 2, 2>>> P_basis_jacobians_transformed_;
  std::vector<DomainType> P_basis_grad_times_u_;
  DomainType omega_minus_xi_D_T_basis_ii_;
  std::vector<double> result1_;
  std::vector<double> result2_;
  std::vector<double> u_coeffs_;
  std::vector<double> P_coeffs_;
  std::vector<double> phi_coeffs_;
  const double frac_c1_Pa_;
  const double grad_u_factor_;
  const double grad_u_T_factor_;
}; // struct OfieldPrepareFunctor

// Computes nonlinear part of ofield residual (-c_1/Pa f(d))
struct CellModelSolver::OfieldResidualFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = OfieldResidualFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  OfieldResidualFunctor(const CellModelSolver& cellmodel,
                        const QuadratureStorage& quad_storage,
                        VectorViewType& res1_vec)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , res1_vec_(res1_vec)
    , P_(0.)
    , result_(quad_storage_.P_basis_size_)
    , P_coeffs_(quad_storage_.P_basis_size_)
    , minus_frac_c1_Pa_(-cellmodel_.c_1_ / cellmodel_.Pa_)
  {}

  OfieldResidualFunctor(const ThisType& other) = default;
  OfieldResidualFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new OfieldResidualFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    std::fill(result_.begin(), result_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_P = cellmodel_.global_indices_P_[index];
    const auto P_basis_size = quad_storage_.P_basis_size_;
    for (size_t ii = 0; ii < P_basis_size; ++ii)
      P_coeffs_[ii] = quad_storage_.P_dofs_[global_indices_P[ii]];
    const auto& quadrature = quad_storage_.quadrature_;
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      P_ *= 0.;
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        P_.axpy(P_coeffs_[ii], quad_storage_.P_basis_evaluated_[ll][ii]);
      const auto factor =
          P_.two_norm2() * minus_frac_c1_Pa_ * quadrature[ll].weight() * quad_storage_.integration_factors_[ll];
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        result_[ii] += (quad_storage_.P_basis_evaluated_[ll][ii] * P_) * factor;
    } // ll
    for (size_t ii = 0; ii < P_basis_size; ++ii)
      res1_vec_.add_to_entry(global_indices_P[ii], result_[ii]);
  }

  const CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  VectorViewType& res1_vec_;
  DomainType P_;
  std::vector<double> result_;
  std::vector<double> P_coeffs_;
  const double minus_frac_c1_Pa_;
}; // struct OfieldResidualFunctor

// Computes nonlinear parts of pfield residual (1/(Be eps^2) f(phi, mu) and 1/eps g(phi))
struct CellModelSolver::PfieldResidualFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = PfieldResidualFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  PfieldResidualFunctor(const CellModelSolver& cellmodel,
                        const QuadratureStorage& quad_storage,
                        VectorViewType& res0_vec,
                        VectorViewType& res1_vec,
                        VectorViewType& res2_vec)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , res0_vec_(res0_vec)
    , res1_vec_(res1_vec)
    , res2_vec_(res2_vec)
    , result1_(quad_storage_.phi_basis_size_)
    , result2_(quad_storage_.phi_basis_size_)
    , phi_coeffs_(quad_storage_.phi_basis_size_)
    , mu_coeffs_(quad_storage_.phi_basis_size_)
    , dt_gamma_div_eps_Ca_(cellmodel_.dt_ * cellmodel_.gamma_ / (cellmodel_.epsilon_ * cellmodel_.Ca_))
    , eps_inv_(1. / cellmodel_.epsilon_)
    , inv_Be_eps2_(1. / (cellmodel_.Be_ * std::pow(cellmodel_.epsilon_, 2)))
  {}

  PfieldResidualFunctor(const ThisType& other) = default;
  PfieldResidualFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new PfieldResidualFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    double phi(0.), mu(0.);
    std::fill(result1_.begin(), result1_.end(), 0.);
    std::fill(result2_.begin(), result2_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const auto phi_basis_size = quad_storage_.phi_basis_size_;
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      phi_coeffs_[ii] = quad_storage_.phi_dofs_[global_indices_phi[ii]];
      mu_coeffs_[ii] = quad_storage_.mu_dofs_[global_indices_phi[ii]];
    }
    const auto& quadrature = quad_storage_.quadrature_;
    const double prefactor2 = cellmodel_.num_pfield_variables_ > 1 ? eps_inv_ : dt_gamma_div_eps_Ca_;
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      phi = XT::Common::transform_reduce(
          phi_coeffs_.begin(), phi_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
      const auto factor = quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      if (cellmodel_.num_pfield_variables_ > 1) {
        mu = XT::Common::transform_reduce(
            mu_coeffs_.begin(), mu_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
        const auto factor1 = inv_Be_eps2_ * (3. * phi * phi - 1) * mu * factor;
        for (size_t ii = 0; ii < phi_basis_size; ++ii)
          result1_[ii] += quad_storage_.phi_basis_evaluated_[ll][ii] * factor1;
      }
      const auto factor2 = prefactor2 * (phi * phi - 1) * phi * factor;
      for (size_t ii = 0; ii < phi_basis_size; ++ii)
        result2_[ii] += quad_storage_.phi_basis_evaluated_[ll][ii] * factor2;
    } // ll
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      const auto index_ii = global_indices_phi[ii];
      if (cellmodel_.num_pfield_variables_ == 1) {
        res0_vec_.add_to_entry(index_ii, result2_[ii]);
      } else {
        res1_vec_.add_to_entry(index_ii, result1_[ii]);
        res2_vec_.add_to_entry(index_ii, result2_[ii]);
      }
    } // ii
  }

  const CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  VectorViewType& res0_vec_;
  VectorViewType& res1_vec_;
  VectorViewType& res2_vec_;
  std::vector<double> result1_;
  std::vector<double> result2_;
  std::vector<double> phi_coeffs_;
  std::vector<double> mu_coeffs_;
  const double dt_gamma_div_eps_Ca_;
  const double eps_inv_;
  const double inv_Be_eps2_;
}; // struct PfieldResidualFunctor

// Computes nonlinear parts of pfield jacobian (1/(Be eps^2) Dphi_f(phi, mu) and Dmu_f(phi) = Dphi_g(phi))
struct CellModelSolver::PfieldNonlinearJacobianFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = PfieldNonlinearJacobianFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  PfieldNonlinearJacobianFunctor(CellModelSolver& cellmodel, const QuadratureStorage& quad_storage)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , result1_(quad_storage_.phi_basis_size_ * quad_storage_.phi_basis_size_)
    , result2_(quad_storage_.phi_basis_size_ * quad_storage_.phi_basis_size_)
    , phi_coeffs_(quad_storage_.phi_basis_size_)
    , mu_coeffs_(quad_storage_.phi_basis_size_)
    , six_inv_Be_eps2_(6. / (cellmodel_.Be_ * std::pow(cellmodel_.epsilon_, 2)))
  {}

  PfieldNonlinearJacobianFunctor(const ThisType& other) = default;
  PfieldNonlinearJacobianFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new PfieldNonlinearJacobianFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    double phi(0.), mu(0.);
    std::fill(result1_.begin(), result1_.end(), 0.);
    std::fill(result2_.begin(), result2_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const auto phi_basis_size = quad_storage_.phi_basis_size_;
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      phi_coeffs_[ii] = quad_storage_.phi_dofs_[global_indices_phi[ii]];
      if (cellmodel_.num_pfield_variables_ > 1)
        mu_coeffs_[ii] = quad_storage_.mu_dofs_[global_indices_phi[ii]];
    }
    const auto& quadrature = quad_storage_.quadrature_;
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      phi = XT::Common::transform_reduce(
          phi_coeffs_.begin(), phi_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
      mu = XT::Common::transform_reduce(
          mu_coeffs_.begin(), mu_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);

      const auto factor = quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      const auto factor1 = (3 * phi * phi - 1) * factor;
      const auto factor2 = six_inv_Be_eps2_ * phi * mu * factor;
      for (size_t ii = 0; ii < phi_basis_size; ++ii) {
        const auto offset_ii = ii * phi_basis_size;
        const auto factor1_ii = quad_storage_.phi_basis_evaluated_[ll][ii] * factor1;
        for (size_t jj = 0; jj < phi_basis_size; ++jj)
          result1_[offset_ii + jj] += quad_storage_.phi_basis_evaluated_[ll][jj] * factor1_ii;
        if (cellmodel_.num_pfield_variables_ > 1) {
          const auto factor2_ii = quad_storage_.phi_basis_evaluated_[ll][ii] * factor2;
          for (size_t jj = 0; jj < phi_basis_size; ++jj)
            result2_[offset_ii + jj] += quad_storage_.phi_basis_evaluated_[ll][jj] * factor2_ii;
        }
      } // ii
    } // ll
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      const auto offset_ii = ii * phi_basis_size;
      const auto index_ii = global_indices_phi[ii];
      for (size_t jj = 0; jj < phi_basis_size; ++jj) {
        const auto index_jj = global_indices_phi[jj];
        cellmodel_.Dmu_f_pfield_.add_to_entry(index_ii, index_jj, result1_[offset_ii + jj]);
        if (cellmodel_.num_pfield_variables_ > 1)
          cellmodel_.Dphi_f_pfield_incl_coeffs_and_sign_.add_to_entry(index_ii, index_jj, result2_[offset_ii + jj]);
      } // jj
    } // ii
  }

  CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  std::vector<double> result1_;
  std::vector<double> result2_;
  std::vector<double> phi_coeffs_;
  std::vector<double> mu_coeffs_;
  const double six_inv_Be_eps2_;
}; // struct PfieldNonlinearJacobianFunctor

// Computes pfield rhs (-c1/(2 Pa) a(d)
struct CellModelSolver::PfieldRhsFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = PfieldRhsFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  PfieldRhsFunctor(const CellModelSolver& cellmodel,
                   const QuadratureStorage& quad_storage,
                   VectorViewType& pfield_r_vector)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , pfield_r_vector_(pfield_r_vector)
    , P_(0.)
    , grad_P_(0.)
    , P_basis_jacobians_transformed_(quad_storage_.P_basis_jacobians_)
    , result_(quad_storage_.phi_basis_size_ * quad_storage_.phi_basis_size_)
    , P_coeffs_(quad_storage_.P_basis_size_)
    , prefactor_((cellmodel_.num_pfield_variables_ == 1 ? cellmodel_.dt_ * cellmodel_.gamma_ : -1.) * cellmodel.c_1_
                 / (2. * cellmodel_.Pa_))
  {}

  PfieldRhsFunctor(const ThisType& other) = default;
  PfieldRhsFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new PfieldRhsFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    std::fill(result_.begin(), result_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_P = cellmodel_.global_indices_P_[index];
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const size_t P_basis_size = quad_storage_.P_basis_size_;
    const size_t phi_basis_size = quad_storage_.phi_basis_size_;
    const auto& quadrature = quad_storage_.quadrature_;
    for (size_t ii = 0; ii < P_basis_size; ++ii)
      P_coeffs_[ii] = quad_storage_.P_dofs_[global_indices_P[ii]];
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      P_ *= 0.;
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        P_.axpy(P_coeffs_[ii], quad_storage_.P_basis_evaluated_[ll][ii]);
      const double factor =
          P_.two_norm2() * prefactor_ * quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      for (size_t ii = 0; ii < phi_basis_size; ++ii)
        result_[ii] += quad_storage_.phi_basis_evaluated_[ll][ii] * factor;
    } // ll
    for (size_t ii = 0; ii < phi_basis_size; ++ii)
      pfield_r_vector_.add_to_entry(global_indices_phi[ii], result_[ii]);
  }

  const CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  VectorViewType& pfield_r_vector_;
  DomainType P_;
  FieldMatrix<double, 2, 2> grad_P_;
  std::vector<std::vector<FieldMatrix<double, 2, 2>>> P_basis_jacobians_transformed_;
  std::vector<double> result_;
  std::vector<double> P_coeffs_;
  const double prefactor_;
  FieldMatrix<double, 2, 2> J_inv_T_;
}; // struct PfieldRhsFunctor

// Computes pfield B(u)
struct CellModelSolver::PfieldBFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = PfieldBFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  PfieldBFunctor(CellModelSolver& cellmodel, const QuadratureStorage& quad_storage)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , u_(0.)
    , J_inv_T_(0.)
    , phi_basis_jacobians_transformed_(quad_storage_.phi_basis_jacobians_)
    , result_(quad_storage_.phi_basis_size_ * quad_storage_.phi_basis_size_)
    , u_coeffs_(quad_storage_.u_basis_size_)
    , phi_coeffs_(quad_storage_.phi_basis_size_)
  {}

  PfieldBFunctor(const ThisType& other) = default;
  PfieldBFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new PfieldBFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const auto& global_indices_u = cellmodel_.global_indices_u_[index];
    const size_t u_basis_size = quad_storage_.u_basis_size_;
    const size_t phi_basis_size = quad_storage_.phi_basis_size_;
    const auto& quadrature = quad_storage_.quadrature_;
    std::fill(result_.begin(), result_.end(), 0.);
    for (size_t ii = 0; ii < u_basis_size; ++ii)
      u_coeffs_[ii] = quad_storage_.u_dofs_[global_indices_u[ii]];
    for (size_t ii = 0; ii < phi_basis_size; ++ii)
      phi_coeffs_[ii] = quad_storage_.phi_dofs_[global_indices_phi[ii]];
    FieldVector<double, 1> basis_grad_ii_times_u;
    const auto& geometry = element.geometry();
    const bool geometry_affine = geometry.affine();
    J_inv_T_ = geometry.jacobianInverseTransposed(geometry.center());
    if (quad_storage_.phi_affine_ && geometry_affine) {
      for (size_t ii = 0; ii < phi_basis_size; ++ii)
        J_inv_T_.mv(quad_storage_.phi_basis_jacobians_[0][ii][0], phi_basis_jacobians_transformed_[0][ii][0]);
    }
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto factor = quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      if (!geometry_affine)
        J_inv_T_ = geometry.jacobianInverseTransposed(quadrature[ll].position());
      if (!quad_storage_.phi_affine_ || !geometry_affine)
        for (size_t ii = 0; ii < phi_basis_size; ++ii)
          J_inv_T_.mv(quad_storage_.phi_basis_jacobians_[ll][ii][0], phi_basis_jacobians_transformed_[ll][ii][0]);

      u_ *= 0.;
      for (size_t ii = 0; ii < u_basis_size; ++ii)
        u_.axpy(u_coeffs_[ii], quad_storage_.u_basis_evaluated_[ll][ii]);
      const auto& phi_basis_grads = quad_storage_.phi_affine_ && geometry_affine ? phi_basis_jacobians_transformed_[0]
                                                                                 : phi_basis_jacobians_transformed_[ll];
      for (size_t ii = 0; ii < phi_basis_size; ++ii) {
        const auto offset_ii = ii * phi_basis_size;
        phi_basis_grads[ii].mv(u_, basis_grad_ii_times_u);
        const double factor_ii = basis_grad_ii_times_u * factor;
        for (size_t jj = 0; jj < phi_basis_size; ++jj)
          result_[offset_ii + jj] += quad_storage_.phi_basis_evaluated_[ll][jj] * factor_ii;
      } // ii
    } // ll
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      const auto offset_ii = ii * phi_basis_size;
      const auto index_ii = global_indices_phi[ii];
      for (size_t jj = 0; jj < phi_basis_size; ++jj)
        cellmodel_.B_pfield_.add_to_entry(index_ii, global_indices_phi[jj], result_[offset_ii + jj]);
    } // ii
  } // ... apply_local(...)

  CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  DomainType u_;
  FieldMatrix<double, 2, 2> J_inv_T_;
  std::vector<std::vector<FieldMatrix<double, 1, d>>> phi_basis_jacobians_transformed_;
  std::vector<double> result_;
  std::vector<double> u_coeffs_;
  std::vector<double> phi_coeffs_;
}; // struct PfieldBFunctor

// Computes Stokes rhs
struct CellModelSolver::StokesRhsFunctor : public XT::Grid::ElementFunctor<PGV>
{
  using ThisType = StokesRhsFunctor;
  using BaseType = typename XT::Grid::ElementFunctor<PGV>;

  StokesRhsFunctor(CellModelSolver& cellmodel, const QuadratureStorage& quad_storage)
    : cellmodel_(cellmodel)
    , quad_storage_(quad_storage)
    , P_(0.)
    , Pnat_(0.)
    , grad_P_(0.)
    , grad_phi_(0.)
    , u_basis_jacobians_transformed_(quad_storage_.u_basis_jacobians_)
    , P_basis_jacobians_transformed_(quad_storage_.P_basis_jacobians_)
    , phi_basis_jacobians_transformed_(quad_storage_.phi_basis_jacobians_)
    , result_(quad_storage_.u_basis_size_)
    , P_coeffs_(quad_storage_.P_basis_size_)
    , Pnat_coeffs_(quad_storage_.P_basis_size_)
    , phi_coeffs_(quad_storage_.phi_basis_size_)
    , phinat_coeffs_(quad_storage_.phi_basis_size_)
    , xi_p_1_(0.5 * (cellmodel_.xi_ + 1))
    , xi_m_1_(0.5 * (cellmodel_.xi_ - 1))
    , eps_div_Ca_(cellmodel_.epsilon_ / cellmodel_.Ca_)
    , Pa_inv_(1. / cellmodel_.Pa_)
  {}

  StokesRhsFunctor(const ThisType& other) = default;
  StokesRhsFunctor(ThisType&&) = default;

  XT::Grid::ElementFunctor<PGV>* copy() final
  {
    return new StokesRhsFunctor(*this);
  }

  void apply_local(const E& element) final
  {
    // TODO: replace 0 by cell index
    std::fill(result_.begin(), result_.end(), 0.);
    const auto index = cellmodel_.grid_view_.indexSet().index(element);
    const auto& global_indices_P = cellmodel_.global_indices_P_[index];
    const auto& global_indices_phi = cellmodel_.global_indices_phi_[index];
    const auto& global_indices_u = cellmodel_.global_indices_u_[index];
    const size_t u_basis_size = quad_storage_.u_basis_size_;
    const size_t P_basis_size = quad_storage_.P_basis_size_;
    const size_t phi_basis_size = quad_storage_.phi_basis_size_;
    const auto& quadrature = quad_storage_.quadrature_;
    for (size_t ii = 0; ii < P_basis_size; ++ii) {
      P_coeffs_[ii] = quad_storage_.P_dofs_[global_indices_P[ii]];
      Pnat_coeffs_[ii] = quad_storage_.Pnat_dofs_[global_indices_P[ii]];
    }
    for (size_t ii = 0; ii < phi_basis_size; ++ii) {
      phi_coeffs_[ii] = quad_storage_.phi_dofs_[global_indices_phi[ii]];
      phinat_coeffs_[ii] = quad_storage_.phinat_dofs_[global_indices_phi[ii]];
    }
    const auto& geometry = element.geometry();
    const bool geometry_affine = geometry.affine();
    J_inv_T_ = geometry.jacobianInverseTransposed(geometry.center());
    if (quad_storage_.P_affine_ && geometry_affine) {
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        for (size_t rr = 0; rr < d; ++rr)
          J_inv_T_.mv(quad_storage_.P_basis_jacobians_[0][ii][rr], P_basis_jacobians_transformed_[0][ii][rr]);
      grad_P_ *= 0.;
      for (size_t ii = 0; ii < P_basis_size; ++ii)
        grad_P_.axpy(P_coeffs_[ii], P_basis_jacobians_transformed_[0][ii]);
    }
    if (quad_storage_.phi_affine_ && geometry_affine) {
      for (size_t ii = 0; ii < phi_basis_size; ++ii)
        J_inv_T_.mv(quad_storage_.phi_basis_jacobians_[0][ii][0], phi_basis_jacobians_transformed_[0][ii][0]);
      grad_phi_ *= 0.;
      for (size_t ii = 0; ii < phi_basis_size; ++ii)
        grad_phi_.axpy(phi_coeffs_[ii], phi_basis_jacobians_transformed_[0][ii]);
    }
    double phi, phinat;
    for (size_t ll = 0; ll < quadrature.size(); ++ll) {
      const auto factor = quad_storage_.integration_factors_[ll] * quadrature[ll].weight();
      if (!geometry_affine)
        J_inv_T_ = geometry.jacobianInverseTransposed(quadrature[ll].position());
      for (size_t ii = 0; ii < u_basis_size; ++ii)
        for (size_t rr = 0; rr < d; ++rr)
          J_inv_T_.mv(quad_storage_.u_basis_jacobians_[ll][ii][rr], u_basis_jacobians_transformed_[ll][ii][rr]);
      if (!quad_storage_.P_affine_ || !geometry_affine)
        for (size_t ii = 0; ii < P_basis_size; ++ii)
          for (size_t rr = 0; rr < d; ++rr)
            J_inv_T_.mv(quad_storage_.P_basis_jacobians_[ll][ii][rr], P_basis_jacobians_transformed_[ll][ii][rr]);
      if (!quad_storage_.phi_affine_ || !geometry_affine)
        for (size_t ii = 0; ii < phi_basis_size; ++ii)
          J_inv_T_.mv(quad_storage_.phi_basis_jacobians_[ll][ii][0], phi_basis_jacobians_transformed_[ll][ii][0]);

      // for (size_t kk = 0; kk < num_cells_; ++kk) {
      P_ *= 0.;
      Pnat_ *= 0.;
      for (size_t ii = 0; ii < P_basis_size; ++ii) {
        P_.axpy(P_coeffs_[ii], quad_storage_.P_basis_evaluated_[ll][ii]);
        Pnat_.axpy(Pnat_coeffs_[ii], quad_storage_.P_basis_evaluated_[ll][ii]);
      }
      if (!quad_storage_.P_affine_ || !geometry_affine) {
        grad_P_ *= 0.;
        for (size_t ii = 0; ii < P_basis_size; ++ii)
          grad_P_.axpy(P_coeffs_[ii], P_basis_jacobians_transformed_[ll][ii]);
      }
      phi = XT::Common::transform_reduce(
          phi_coeffs_.begin(), phi_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
      phinat = XT::Common::transform_reduce(
          phinat_coeffs_.begin(), phinat_coeffs_.end(), quad_storage_.phi_basis_evaluated_[ll].begin(), 0.);
      if (!quad_storage_.phi_affine_ || !geometry_affine) {
        grad_phi_ *= 0.;
        for (size_t ii = 0; ii < phi_basis_size; ++ii)
          grad_phi_.axpy(phi_coeffs_[ii], phi_basis_jacobians_transformed_[ll][ii]);
      }
      const auto phi_tilde = (phi + 1.) / 2.;

      // evaluate rhs terms
      auto grad_P_T_Pnat_plus_phinat_grad_phi = Pnat_;
      grad_P_.mtv(Pnat_, grad_P_T_Pnat_plus_phinat_grad_phi);
      grad_P_T_Pnat_plus_phinat_grad_phi.axpy(phinat, grad_phi_[0]);
      const auto Fa_inv_phi_tilde = cellmodel_.Fa_inv_ * phi_tilde;
      for (size_t ii = 0; ii < u_basis_size; ++ii) {
        const auto& u_basis_grad_ii = u_basis_jacobians_transformed_[ll][ii];
        const auto& u_basis_val_ii = quad_storage_.u_basis_evaluated_[ll][ii];
        for (size_t mm = 0; mm < d; ++mm) {
          const auto factor1_mm = (Fa_inv_phi_tilde * P_[mm] + xi_p_1_ * Pnat_[mm]) * factor;
          const auto factor2_mm = xi_m_1_ * P_[mm] * factor;
          if (cellmodel_.num_pfield_variables_ > 1)
            result_[ii] += grad_P_T_Pnat_plus_phinat_grad_phi[mm] * u_basis_val_ii[mm] * factor;
          result_[ii] += -(P_ * u_basis_grad_ii[mm]) * factor1_mm - (Pnat_ * u_basis_grad_ii[mm]) * factor2_mm;
          if (cellmodel_.num_pfield_variables_ == 1) {
            const auto factor3_mm = (eps_div_Ca_ * grad_phi_[0][mm]) * factor;
            result_[ii] += (grad_phi_[0] * u_basis_grad_ii[mm]) * factor3_mm;
            for (size_t nn = 0; nn < d; ++nn) {
              double grad_P_T_grad_P_mn = 0.;
              for (size_t kk = 0; kk < d; ++kk)
                grad_P_T_grad_P_mn += grad_P_[mm][kk] * grad_P_[kk][nn];
              result_[ii] += grad_P_T_grad_P_mn * u_basis_grad_ii[mm][nn] * Pa_inv_ * factor;
            } // nn
          } // if (num_pfield_variables_ == 1)
        } // mm
      } // ii
    } // ll
    for (size_t ii = 0; ii < u_basis_size; ++ii)
      cellmodel_.stokes_f_vector_.add_to_entry(global_indices_u[ii], result_[ii]);
  }

  CellModelSolver& cellmodel_;
  const QuadratureStorage& quad_storage_;
  DomainType P_;
  DomainType Pnat_;
  FieldMatrix<double, 2, 2> grad_P_;
  FieldMatrix<double, 1, 2> grad_phi_;
  std::vector<std::vector<FieldMatrix<double, 2, 2>>> u_basis_jacobians_transformed_;
  std::vector<std::vector<FieldMatrix<double, 2, 2>>> P_basis_jacobians_transformed_;
  std::vector<std::vector<FieldMatrix<double, 1, 2>>> phi_basis_jacobians_transformed_;
  std::vector<double> result_;
  std::vector<double> P_coeffs_;
  std::vector<double> Pnat_coeffs_;
  std::vector<double> phi_coeffs_;
  std::vector<double> phinat_coeffs_;
  const double xi_p_1_;
  const double xi_m_1_;
  const double eps_div_Ca_;
  const double Pa_inv_;
  FieldMatrix<double, 2, 2> J_inv_T_;
}; // struct StokesRhsFunctor


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
 *Ignored if write = false. filename: Prefix for .vtu and .txt files. Ignored if write = false. subsampling: Whether
 *to use subsampling for visualization. Ignored if write = false.
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

void CellModelSolver::solve_without_storing(const bool write,
                                            const double write_step,
                                            const std::string filename,
                                            const bool subsampling)
{
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
      pfield_vectors_[kk] = apply_inverse_pfield_operator(pfield_vectors_[kk], kk);
      std::cout << "Pfield " << kk << " done" << std::endl;
      prepare_ofield_operator(kk);
      ofield_vectors_[kk] = apply_inverse_ofield_operator(ofield_vectors_[kk], kk);
      std::cout << "Ofield " << kk << " done" << std::endl;
    }

    // stokes system
    prepare_stokes_operator();
    stokes_vector_ = apply_inverse_stokes_operator();
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
  VectorType ret(num_pfield_variables_ * size_phi_);
  ConstVectorViewType phi_view(vec, 0, size_phi_);
  VectorViewType phi_ret_view(ret, 0, size_phi_);
  M_pfield_.mv(phi_view, phi_ret_view);
  if (num_pfield_variables_ > 1) {
    ConstVectorViewType phinat_view(vec, size_phi_, 2 * size_phi_);
    ConstVectorViewType mu_view(vec, 2 * size_phi_, 3 * size_phi_);
    VectorViewType phinat_ret_view(ret, size_phi_, 2 * size_phi_);
    VectorViewType mu_ret_view(ret, 2 * size_phi_, 3 * size_phi_);
    M_pfield_.mv(phinat_view, phinat_ret_view);
    M_pfield_.mv(mu_view, mu_ret_view);
  }
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
  // M_p_stokes_.mv(p_view, p_ret_view);
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
  const auto phi_func = make_discrete_function(phi_space_, phi_vec, "phi");
  XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phi_func);
  if (num_pfield_variables_ > 1) {
    const ConstVectorViewType phinat_vec(vec, size_phi_, 2 * size_phi_);
    const ConstVectorViewType mu_vec(vec, 2 * size_phi_, 3 * size_phi_);
    const auto phinat_func = make_discrete_function(phi_space_, phinat_vec, "phinat");
    const auto mu_func = make_discrete_function(phi_space_, mu_vec, "mu");
    XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phinat_func);
    XT::Functions::internal::add_to_vtkwriter(*vtk_writer, mu_func);
  }
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
                                const bool txt,
                                const bool timings) const
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
      if (num_pfield_variables_ > 1) {
        XT::Functions::internal::add_to_vtkwriter(*vtk_writer, phinat_[kk]);
        XT::Functions::internal::add_to_vtkwriter(*vtk_writer, mu_[kk]);
      }
      // XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, *phi_funcs[kk]);
      XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, phi_[kk]);
      if (num_pfield_variables_ > 1) {
        XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, phinat_[kk]);
        XT::Functions::internal::add_gradient_to_vtkwriter(*vtk_writer, mu_[kk]);
      }
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
      if (num_pfield_variables_ > 1) {
        write_to_textfile(phinat_[kk], prefix, step, t);
        write_to_textfile(mu_[kk], prefix, step, t);
      }
    } // kk
  } // if (txt)
  if (timings) {
    const auto& pfield_solve_times = pfield_solver_statistics_.solve_time_;
    const auto& ofield_solve_times = ofield_solver_statistics_.solve_time_;
    const auto& stokes_solve_times = stokes_solver_statistics_.solve_time_;
    const auto& pfield_iterations = pfield_solver_statistics_.iterations_;
    const auto& ofield_iterations = ofield_solver_statistics_.iterations_;
    const auto& stokes_iterations = stokes_solver_statistics_.iterations_;
    assert(pfield_solve_times.size() == pfield_iterations.size());
    assert(ofield_solve_times.size() == ofield_iterations.size());
    assert(stokes_solve_times.size() == stokes_iterations.size());
    const double pfield_total_solve_time =
        std::accumulate(pfield_solve_times.begin(), pfield_solve_times.end(), std::chrono::duration<double>(0.))
            .count();
    const double ofield_total_solve_time =
        std::accumulate(ofield_solve_times.begin(), ofield_solve_times.end(), std::chrono::duration<double>(0.))
            .count();
    const double stokes_total_solve_time =
        std::accumulate(stokes_solve_times.begin(), stokes_solve_times.end(), std::chrono::duration<double>(0.))
            .count();
    const double pfield_average_time_per_solve =
        pfield_solve_times.size() == 0 ? 0. : pfield_total_solve_time / pfield_solve_times.size();
    const double ofield_average_time_per_solve =
        ofield_solve_times.size() == 0 ? 0. : ofield_total_solve_time / ofield_solve_times.size();
    const double stokes_average_time_per_solve =
        stokes_solve_times.size() == 0 ? 0. : stokes_total_solve_time / stokes_solve_times.size();
    const double pfield_average_num_iterations =
        pfield_solve_times.size() == 0
            ? 0.
            : std::accumulate(pfield_iterations.begin(), pfield_iterations.end(), 0) / double(pfield_iterations.size());
    const double ofield_average_num_iterations =
        ofield_solve_times.size() == 0
            ? 0.
            : std::accumulate(ofield_iterations.begin(), ofield_iterations.end(), 0) / double(ofield_iterations.size());
    const double stokes_average_num_iterations =
        stokes_solve_times.size() == 0
            ? 0.
            : std::accumulate(stokes_iterations.begin(), stokes_iterations.end(), 0) / double(stokes_iterations.size());
    std::string timings_filename = prefix + "_timings.txt";
    // check if file exists
    bool file_exists;
    {
      std::ifstream timings_file(timings_filename);
      file_exists = timings_file.good();
    }
    std::ofstream timings_file(timings_filename, std::ios_base::app);
    if (!file_exists) {
      timings_file << "pf_solves pf_average_its pf_setup_time pf_av_solve_time | "
                   << "of_solves of_average_its of_setup_time of_av_solve_time | "
                   << "st_solves st_average_its st_setup_time st_av_solve_time" << std::endl;
    }
    timings_file << pfield_iterations.size() << " " << pfield_average_num_iterations << " " << std::fixed
                 << std::setprecision(4) << pfield_solver_statistics_.setup_time_.count() << " "
                 << pfield_average_time_per_solve << " | " << ofield_iterations.size() << " "
                 << ofield_average_num_iterations << " " << ofield_solver_statistics_.setup_time_.count() << " "
                 << ofield_average_time_per_solve << " | " << stokes_iterations.size() << " "
                 << stokes_average_num_iterations << " " << stokes_solver_statistics_.setup_time_.count() << " "
                 << stokes_average_time_per_solve << std::endl;
    timings_file.close();
  }
} // void visualize(...)

//******************************************************************************************************************
//*******************************  Methods to get and set variable values   ****************************************
//******************************************************************************************************************

// Sets stokes vector to stokes_vec
void CellModelSolver::set_stokes_vec(const VectorType& stokes_vec)
{
#ifndef NDEBUG
  DUNE_THROW_IF(
      stokes_vec.size() != size_u_ + size_p_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
#endif
  stokes_vector_ = stokes_vec;
}

// Sets given dofs of stokes vector to values
void CellModelSolver::set_stokes_vec_dofs(const std::vector<R>& values, const std::vector<size_t>& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    stokes_vector_.set_entry(dofs[ii], values[ii]);
}

// Sets given dofs of stokes vector to values
void CellModelSolver::set_stokes_vec_dofs(const pybind11::array_t<double>& values, const pybind11::list& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(static_cast<size_t>(values.size()) != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  const auto& values_data = values.unchecked();
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    stokes_vector_.set_entry(dofs[ii].cast<size_t>(), values_data[ii]);
}

// Sets orientation field vector belonging to cell to pfield_vec
void CellModelSolver::set_ofield_vec(const size_t cell, const VectorType& ofield_vec)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(ofield_vec.size() != 2 * size_P_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
#endif
  ofield_vectors_[cell] = ofield_vec;
}

// Sets given dofs of orientation field vector belonging to cell to values
void CellModelSolver::set_ofield_vec_dofs(const size_t cell,
                                          const std::vector<R>& values,
                                          const std::vector<size_t>& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    ofield_vectors_[cell].set_entry(dofs[ii], values[ii]);
}

// Sets given dofs of orientation field vector belonging to cell to values
void CellModelSolver::set_ofield_vec_dofs(const size_t cell,
                                          const pybind11::array_t<double>& values,
                                          const pybind11::list& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(static_cast<size_t>(values.size()) != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  const auto& values_data = values.unchecked();
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    ofield_vectors_[cell].set_entry(dofs[ii].cast<size_t>(), values_data[ii]);
}

// Sets phasefield vector belonging to cell to pfield_vec
void CellModelSolver::set_pfield_vec(const size_t cell, const VectorType& pfield_vec)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(pfield_vec.size() != 3 * size_phi_, XT::Common::Exceptions::wrong_input_given, "Invalid vector size!");
#endif
  pfield_vectors_[cell] = pfield_vec;
}

// Sets given dofs of phasefield vector belonging to cell to values
void CellModelSolver::set_pfield_vec_dofs(const size_t cell,
                                          const std::vector<R>& values,
                                          const std::vector<size_t>& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(values.size() != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    pfield_vectors_[cell].set_entry(dofs[ii], values[ii]);
}

// Sets given dofs of phasefield vector belonging to cell to values
void CellModelSolver::set_pfield_vec_dofs(const size_t cell,
                                          const pybind11::array_t<double>& values,
                                          const pybind11::list& dofs)
{
#ifndef NDEBUG
  DUNE_THROW_IF(cell >= num_cells_, XT::Common::Exceptions::wrong_input_given, "Invalid cell index");
  DUNE_THROW_IF(static_cast<size_t>(values.size()) != dofs.size(),
                XT::Common::Exceptions::wrong_input_given,
                "Size of values does not match size of dofs");
#endif
  const auto& values_data = values.unchecked();
  for (size_t ii = 0; ii < dofs.size(); ++ii)
    pfield_vectors_[cell].set_entry(dofs[ii].cast<size_t>(), values_data[ii]);
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
  // copy VectorViews to actual vectors to improve performance
  u_tmp_.dofs().vector() = u_.dofs().vector();
  for (size_t kk = 0; kk < num_cells_; kk++) {
    phi_tmp_[kk].dofs().vector() = phi_[kk].dofs().vector();
    if (num_pfield_variables_ > 1)
      phinat_tmp_[kk].dofs().vector() = phinat_[kk].dofs().vector();
    P_tmp_[kk].dofs().vector() = P_[kk].dofs().vector();
    Pnat_tmp_[kk].dofs().vector() = Pnat_[kk].dofs().vector();
  }
  // assemble
  assemble_stokes_rhs(restricted);
  stokes_jac_linear_op_.prepare(0, restricted);
}

void CellModelSolver::prepare_ofield_operator(const size_t cell, const bool restricted)
{
  // copy VectorViews to actual vectors to improve performance
  u_tmp_.dofs().vector() = u_.dofs().vector();
  P_tmp_[cell].dofs().vector() = P_[cell].dofs().vector();
  phi_tmp_[cell].dofs().vector() = phi_[cell].dofs().vector();
  // assemble
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
    if (num_pfield_variables_ > 1)
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
    p_deim_source_dofs_begin_ =
        std::lower_bound(source_dofs[2].begin(), source_dofs[2].end(), size_u_) - source_dofs[2].begin();
    // extract relevant cols of B^T to be able to apply B rapidly
    B_stokes_restricted_ = DenseMatrixType(p_range_dofs.size(), p_deim_source_dofs_begin_, 0., 0);
    for (size_t ii = 0; ii < p_range_dofs.size(); ++ii)
      for (size_t jj = 0; jj < p_deim_source_dofs_begin_; ++jj)
        B_stokes_restricted_.set_entry(ii, jj, BT_stokes_.get_entry(source_dofs[2][jj], p_range_dofs[ii]));
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
    if (num_pfield_variables_ != 3)
      DUNE_THROW(Dune::NotImplemented, "Restricted operator has not yet been adapted!");
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
    get_deim_source_dofs(source_dofs[0],
                         pfield_solver_.system_matrix_pattern(pfield_submatrix_pattern_, num_pfield_variables_),
                         unique_range_dofs);

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
  auto& residual = stokes_tmp_vec2_;
  // copy values to high-dimensional vector
  stokes_jac_linear_op_.apply(y, residual);
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
  }
  return residual;
}

void CellModelSolver::assemble_nonlinear_part_of_ofield_residual(VectorType& residual,
                                                                 const size_t cell,
                                                                 const bool restricted)
{
  // nonlinear part
  VectorViewType res1_vec(residual, size_P_, 2 * size_P_);
  OfieldResidualFunctor functor(*this, *ofield_residual_and_nonlinear_S10_quad_, res1_vec);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(ofield_deim_entities_[cell]);
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
    mv(Dd_f_ofield_incl_coeffs_and_sign_, source_P, P_tmp_vec_, Pnat_range_dofs);
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
  }
  return residual;
}

void CellModelSolver::update_ofield_parameters(const double Pa, const size_t /*cell*/, const bool restricted)
{
  Pa_ = Pa;
  if (!restricted)
    ofield_solver_.set_params(dt_, kappa_, Pa_);
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
    const ConstVectorViewType source_phi(source, 0, size_phi_);
    // nonlinear_part
    auto& tmp_vec = phi_tmp_vec_;
    const auto mv = mv_func<ConstVectorViewType, VectorType>(restricted);
    const auto vector_axpy = vector_axpy_func<VectorViewType, VectorType>(restricted);
    const auto scal = scal_func<VectorType>(restricted);
    const auto add = add_func<VectorViewType, VectorType>(restricted);
    if (num_pfield_variables_ == 1) {
      const auto& phi_range_dofs = phi_deim_range_dofs_[cell];
      VectorViewType range_phi(residual, 0, size_phi_);
      mv(Dmu_f_pfield_, source_phi, tmp_vec, phi_range_dofs);
      vector_axpy(range_phi, dt_ * gamma_ / (epsilon_ * Ca_), tmp_vec, phi_range_dofs);
    } else {
      const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
      const auto& mu_range_dofs = mu_deim_range_dofs_[cell];
      VectorViewType range_phinat(residual, size_phi_, 2 * size_phi_);
      VectorViewType range_mu(residual, 2 * size_phi_, 3 * size_phi_);
      const ConstVectorViewType source_mu(source, 2 * size_phi_, 3 * size_phi_);
      // apply missing parts of S_12
      mv(Dmu_f_pfield_, source_mu, tmp_vec, phinat_range_dofs);
      vector_axpy(range_phinat, 1. / (Be_ * std::pow(epsilon_, 2)), tmp_vec, phinat_range_dofs);
      // apply missing parts of S_20
      mv(Dmu_f_pfield_, source_phi, tmp_vec, mu_range_dofs);
      vector_axpy(range_mu, 1. / epsilon_, tmp_vec, mu_range_dofs);
      // apply S_10
      mv(Dphi_f_pfield_incl_coeffs_and_sign_, source_phi, tmp_vec, phinat_range_dofs);
      add(range_phinat, tmp_vec, phinat_range_dofs);
    }
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
  }
  return residual;
}

void CellModelSolver::update_pfield_parameters(
    const double Be, const double Ca, const double Pa, const size_t /*cell*/, const bool restricted)
{
  Be_ = Be;
  Ca_ = Ca;
  Pa_ = Pa;
  if (!restricted)
    pfield_solver_.set_params(gamma_, epsilon_, Be_, Ca_);
}

//******************************************************************************************************************
//******************************************* Apply inverse operators **********************************************
//******************************************************************************************************************

// Applies inverse stokes operator (solves F(y) = 0)
CellModelSolver::VectorType CellModelSolver::apply_inverse_stokes_operator()
{
  return apply_inverse_stokes_helper(stokes_rhs_vector_);
}

CellModelSolver::VectorType CellModelSolver::apply_inverse_stokes_helper(const EigenVectorType& rhs)
{
  EigenVectorType ret(size_u_ + size_p_);
  auto t1 = std::chrono::high_resolution_clock::now();
  if (stokes_solver_type_ == StokesSolverType::direct) {
    ret.backend() = stokes_solver_->solve(rhs.backend());
    stokes_solver_statistics_.iterations_.emplace_back(1);
    stokes_solver_statistics_.solve_time_.emplace_back(std::chrono::high_resolution_clock::now() - t1);
    return XT::Common::convert_to<VectorType>(ret);
  }
  ConstEigenVectorViewType rhs_u_view(rhs, 0, size_u_);
  // calculate rhs B A^{-1} f
  auto& p = stokes_p_tmp_vec_;
  auto& B_Ainv_f = stokes_p_tmp_vec2_;
  auto& Ainv_f = stokes_u_tmp_vec_;
  auto& rhs_u = stokes_u_tmp_vec2_;
  rhs_u = rhs_u_view;
  Ainv_f.backend() = stokes_A_solver_->solve(rhs_u.backend());
  BT_stokes_.mtv(Ainv_f, B_Ainv_f);
  // Solve S p = rhs
  InverseOperatorResult res;
  stokes_schur_solver_->apply(p, B_Ainv_f, res);
  // Now solve u = A^{-1}(f - B^T p)
  auto& BT_p = stokes_u_tmp_vec_;
  BT_stokes_.mv(p, BT_p);
  // rhs_u already contains f
  rhs_u -= BT_p;
  auto& u = stokes_u_tmp_vec_;
  u.backend() = stokes_A_solver_->solve(rhs_u.backend());
  for (size_t ii = 0; ii < size_u_; ++ii)
    ret.set_entry(ii, u.get_entry(ii));
  for (size_t ii = 0; ii < size_p_; ++ii)
    ret.set_entry(size_u_ + ii, p.get_entry(ii));
  stokes_solver_statistics_.iterations_.emplace_back(res.iterations);
  stokes_solver_statistics_.solve_time_.emplace_back(std::chrono::high_resolution_clock::now() - t1);
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
    auto t1 = std::chrono::high_resolution_clock::now();
    update = ofield_solver_.apply(residual, cell);
    ofield_solver_statistics_.iterations_.emplace_back(ofield_solver_.wrapper_.statistics_.iterations);
    ofield_solver_statistics_.solve_time_.emplace_back(std::chrono::high_resolution_clock::now() - t1);

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
    auto t1 = std::chrono::high_resolution_clock::now();
    update = pfield_solver_.apply(residual, cell);
    pfield_solver_statistics_.iterations_.emplace_back(pfield_solver_.wrapper_.statistics_.iterations);
    pfield_solver_statistics_.solve_time_.emplace_back(std::chrono::high_resolution_clock::now() - t1);

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
  return apply_inverse_stokes_helper(rhs_eigen);
}

//******************************************************************************************************************
//**************************** Methods to assemble rhs, residuals and jacobians ************************************
//******************************************************************************************************************

// Computes stokes rhs using currently stored values of variables and stores in stokes_rhs_vector_
void CellModelSolver::assemble_stokes_rhs(const bool restricted)
{
  const auto scal = scal_func<EigenVectorViewType>(restricted);
  const auto& u_range_dofs = u_deim_range_dofs_;
  scal(stokes_f_vector_, 0., u_range_dofs);
  StokesRhsFunctor functor(*this, *stokes_rhs_quad_);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(stokes_deim_entities_);
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
}

// Computes phase field rhs using currently stored values of variables and stores in pfield_rhs_vector_
void CellModelSolver::assemble_pfield_rhs(const size_t cell, const bool restricted)
{
  const auto& phi_range_dofs = phi_deim_range_dofs_[cell];
  // calculate r0
  const auto mv = mv_func<VectorType, VectorViewType>(restricted);
  auto& phi_n = phi_tmp_vec2_;
  phi_n = phi_[cell].dofs().vector();
  mv(M_pfield_, phi_n, pfield_r0_vector_, phi_range_dofs);
  XT::Grid::Walker<PGV> walker(grid_view_);
  const auto scal = scal_func<VectorViewType>(restricted);
  if (num_pfield_variables_ == 1) {
    PfieldRhsFunctor functor(*this, *pfield_rhs_quad_, pfield_r0_vector_);
    walker.append(functor);
  } else {
    const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
    scal(pfield_r1_vector_, 0., phinat_range_dofs);
    PfieldRhsFunctor functor(*this, *pfield_rhs_quad_, pfield_r1_vector_);
    walker.append(functor);
  }
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(pfield_deim_entities_[cell]);
}

// assembles linear part of orientation field jacobian and stores in S_ofield_
void CellModelSolver::assemble_ofield_linear_jacobian(const size_t cell, const bool restricted)
{
  const auto& P_range_dofs = P_deim_range_dofs_[cell];
  const auto& Pnat_range_dofs = Pnat_deim_range_dofs_[cell];
  set_mat_to_zero(B_ofield_, restricted, ofield_submatrix_pattern_, P_range_dofs);
  set_mat_to_zero(C_ofield_incl_coeffs_and_sign_, restricted, ofield_submatrix_pattern_, Pnat_range_dofs);
  OfieldPrepareFunctor functor(*this, *ofield_prepare_quad_);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(ofield_deim_entities_[cell]);
  // In the restricted case, we cannot use the solver anyway
  if (ofield_solver_.is_schur_solver() && !restricted) {
    S_schur_ofield_linear_part_.backend() = M_ofield_.backend();
    S_schur_ofield_linear_part_.axpy(dt_, B_ofield_);
    S_schur_ofield_linear_part_.axpy(-dt_ / kappa_, C_ofield_incl_coeffs_and_sign_);
    S_schur_ofield_linear_part_.axpy(dt_ / (kappa_ * Pa_), E_ofield_);
  }
} // ... assemble_ofield_linear_jacobian(...)

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
      ofield_solver_.schur_matrix().axpy(-dt_ / kappa_, Dd_f_ofield_incl_coeffs_and_sign_);
    }
    ofield_solver_.prepare(cell, restricted);
  }
}

void CellModelSolver::assemble_C_ofield_nonlinear_part(const size_t cell, const bool restricted)
{
  set_mat_to_zero(
      Dd_f_ofield_incl_coeffs_and_sign_, restricted, ofield_submatrix_pattern_, Pnat_deim_range_dofs_[cell]);
  OfieldNonlinearS10Functor functor(*this, *ofield_residual_and_nonlinear_S10_quad_);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(ofield_deim_entities_[cell]);
}

// assembles linear part of phase field jacobian
void CellModelSolver::assemble_pfield_linear_jacobian(const size_t cell, const bool restricted)
{
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
  PfieldBFunctor functor(*this, *pfield_B_quad_);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(pfield_deim_entities_[cell]);
}

// assembles nonlinear part of phase field jacobian
void CellModelSolver::assemble_pfield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted)
{
  fill_tmp_pfield(cell, y, restricted);
  // clear matrices
  if (num_pfield_variables_ == 1) {
    set_mat_to_zero(Dmu_f_pfield_, restricted, pfield_submatrix_pattern_, phi_deim_range_dofs_[cell]);
  } else {
    const auto& phinat_range_dofs = phinat_deim_range_dofs_[cell];
    const auto& mu_range_dofs = mu_deim_range_dofs_[cell];
    set_mat_to_zero(Dmu_f_pfield_, restricted, pfield_submatrix_pattern_, mu_range_dofs);
    if (restricted)
      set_mat_to_zero(Dmu_f_pfield_, restricted, pfield_submatrix_pattern_, phinat_range_dofs);
    set_mat_to_zero(Dphi_f_pfield_incl_coeffs_and_sign_, restricted, pfield_submatrix_pattern_, phinat_range_dofs);
  }
  DUNE_THROW_IF(num_cells_ > 1, Dune::NotImplemented, "");
  PfieldNonlinearJacobianFunctor functor(*this, *pfield_residual_and_nonlinear_jac_quad_);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(pfield_deim_entities_[cell]);
}

// assembles nonlinear part of phasefield residual and adds to residual
void CellModelSolver::assemble_nonlinear_part_of_pfield_residual(VectorType& residual,
                                                                 const size_t cell,
                                                                 const bool restricted)
{
  VectorViewType res0_vec(residual, 0, size_phi_);
  VectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
  VectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);
  DUNE_THROW_IF(num_cells_ > 1, Dune::NotImplemented, "");
  PfieldResidualFunctor functor(*this, *pfield_residual_and_nonlinear_jac_quad_, res0_vec, res1_vec, res2_vec);
  XT::Grid::Walker<PGV> walker(grid_view_);
  walker.append(functor);
  if (!restricted)
    walker.walk(partitioning_);
  // walker.walk(use_tbb_);
  else
    walker.walk_range(pfield_deim_entities_[cell]);
} // ... assemble_nonlinear_part_of_pfield_residual(...)

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
  if (testcase == "single_cell" || testcase == "single_cell_dirichlet" || testcase == "channel"
      || testcase == "cell_isolation_experiment")
    return {{0., 0.}};
  DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return FieldVector<R, d>();
}

// get upper right of computational domain from testcase name
XT::Common::FieldVector<CellModelSolver::R, CellModelSolver::d>
CellModelSolver::get_upper_right(const std::string& testcase)
{
  if (testcase == "single_cell" || testcase == "single_cell_dirichlet")
    return {{30., 30.}};
  if (testcase == "cell_isolation_experiment")
    return {{40., 40.}};
  if (testcase == "channel")
    return {{160., 40.}};
  DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return FieldVector<R, d>();
}

// get directions in which domain is periodic from testcase name
std::string CellModelSolver::get_periodic_directions(const std::string& testcase)
{
  if (testcase == "single_cell" || testcase == "channel")
    return "01";
  if (testcase == "single_cell_dirichlet" || testcase == "cell_isolation_experiment")
    return "00";
  DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
  return "";
}

// get number of cells from testcase name
size_t CellModelSolver::get_num_cells(const std::string& testcase)
{
  if (testcase == "single_cell" || testcase == "single_cell_dirichlet" || testcase == "channel"
      || testcase == "cell_isolation_experiment")
    return 1;
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

// Copies P_, Pnat_ (which are VectorViews) to vectors to improve performance in the subsequent operations
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

// Copies phi_, phinat_, mu_ (which are VectorViews) to vectors to improve performance in the subsequent operations
void CellModelSolver::fill_tmp_pfield(const size_t cell, const VectorType& source, const bool restricted) const
{
  if (!restricted) {
    ConstVectorViewType phi_vec(source, 0, size_phi_);
    phi_tmp_[cell].dofs().vector() = phi_vec;
    if (num_pfield_variables_ > 1) {
      ConstVectorViewType mu_vec(source, 2 * size_phi_, 3 * size_phi_);
      mu_tmp_[cell].dofs().vector() = mu_vec;
    }
  } else {
    if (num_pfield_variables_ == 1)
      DUNE_THROW(Dune::NotImplemented, "Restricted operator has not yet been adapted!");
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
