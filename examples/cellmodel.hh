// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_EXAMPLES_CELLMODEL_HH
#define DUNE_GDT_EXAMPLES_CELLMODEL_HH

#include <chrono>
#include <cstdlib>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/common.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/container/matrix-view.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/grid/boundaryinfo.hh>
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
#include <dune/gdt/local/integrands/symmetric_elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/integrands/gradient_value.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/norms.hh>

#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>

#include <Eigen/IterativeLinearSolvers>

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

struct CellModelSolver
{
  // using G = ALU_2D_SIMPLEX_CONFORMING;
  using G = YASP_2D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  using GV = typename G::LeafGridView;
  using PGV = XT::Grid::PeriodicGridView<GV>;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using PI = XT::Grid::extract_intersection_t<PGV>;
  using MatrixType = XT::LA::EigenRowMajorSparseMatrix<double>;
  //   using VectorType = XT::LA::EigenDenseVector<double>;
  using VectorType = XT::LA::CommonDenseVector<double>;
  using EigenVectorType = XT::LA::EigenDenseVector<double>;
  // using MatrixType = XT::LA::IstlRowMajorSparseMatrix<double>;
  // using VectorType = XT::LA::IstlDenseVector<double>;
  using MatrixViewType = XT::LA::MatrixView<MatrixType>;
  using VectorViewType = XT::LA::VectorView<VectorType>;
  using EigenVectorViewType = XT::LA::VectorView<EigenVectorType>;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;
  using ConstEigenVectorViewType = XT::LA::ConstVectorView<EigenVectorType>;
  using R = typename XT::Functions::GenericGridFunction<E, d>::RangeFieldType;
  using DiscreteFunctionType = DiscreteFunction<VectorType, PGV, 1, 1, R>;
  using LocalDiscreteFunctionType = typename DiscreteFunctionType::LocalFunctionType;
  using VectorDiscreteFunctionType = DiscreteFunction<VectorType, PGV, d, 1, R>;
  using VectorLocalDiscreteFunctionType = typename VectorDiscreteFunctionType::LocalFunctionType;
  using ViewDiscreteFunctionType = DiscreteFunction<VectorViewType, PGV, 1, 1, R>;
  using ViewLocalDiscreteFunctionType = typename ViewDiscreteFunctionType::LocalFunctionType;
  using ViewVectorDiscreteFunctionType = DiscreteFunction<VectorViewType, PGV, d, 1, R>;
  using ViewVectorLocalDiscreteFunctionType = typename ViewVectorDiscreteFunctionType::LocalFunctionType;
  using DomainType = typename XT::Functions::GenericGridFunction<E, d>::DomainType;
  using DomainRetType = XT::Common::FieldVector<R, d>;
  using JacobianRetType = XT::Common::FieldMatrix<R, d, d>;
  using ColMajorBackendType = ::Eigen::SparseMatrix<R, ::Eigen::ColMajor>;
  using RowMajorBackendType = typename MatrixType::BackendType;
  using StokesSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using SolverType = ::Eigen::BiCGSTAB<RowMajorBackendType, ::Eigen::IncompleteLUT<R>>;
  using StokesPerThread = XT::Common::PerThreadValue<std::unique_ptr<VectorLocalDiscreteFunctionType>>;
  using PfieldPerThread = XT::Common::PerThreadValue<std::vector<std::unique_ptr<LocalDiscreteFunctionType>>>;
  using OfieldPerThread = XT::Common::PerThreadValue<std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>>>;

  CellModelSolver(const std::string testcase = "single_cell",
                  const double t_end = 1.,
                  const unsigned int num_elements_x = 50,
                  const unsigned int num_elements_y = 50,
                  const double Re = 5e-13,
                  const double Fa = 1.,
                  const double xi = 1.1,
                  const double kappa = 1.65,
                  const double c_1 = 5.,
                  const double Pa = 1,
                  const double beta = 0.,
                  const double gamma = 0.025,
                  const double Be = 0.3,
                  const double Ca = 0.1,
                  const double epsilon = 0.21,
                  const double In = 1.,
                  const bool linearize = false,
                  const double pol_order = 2)
    // do a global refine once, this makes simplicial grids look more symmetric
    : lower_left_(get_lower_left(testcase))
    , upper_right_(get_upper_right(testcase))
    , grid_(XT::Grid::make_cube_grid<G>(lower_left_, upper_right_, {num_elements_x, num_elements_y}, 1))
    , t_end_(t_end)
    , t_(0.)
    , nonperiodic_grid_view_(grid_.leaf_view())
    , grid_view_(nonperiodic_grid_view_, std::bitset<d>(get_periodic_directions(testcase)))
    , u_space_(make_continuous_lagrange_space<d>(grid_view_, pol_order))
    , p_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order - 1))
    , phi_space_(make_continuous_lagrange_space<1>(grid_view_, pol_order))
    , size_u_(u_space_.mapper().size())
    , size_p_(p_space_.mapper().size())
    , size_phi_(phi_space_.mapper().size())
    , num_cells_(get_num_cells(testcase))
    , stokes_vector_(size_u_ + size_p_, 0.)
    , ofield_vector_(2 * size_u_, 0.)
    , pfield_vector_(3 * size_phi_, 0.)
    , u_view_(stokes_vector_, 0, size_u_)
    , p_view_(stokes_vector_, size_u_, size_u_ + size_p_)
    , P_view_(ofield_vector_, 0, size_u_)
    , Pnat_view_(ofield_vector_, size_u_, 2 * size_u_)
    , phi_view_(pfield_vector_, 0, size_phi_)
    , phinat_view_(pfield_vector_, size_phi_, 2 * size_phi_)
    , mu_view_(pfield_vector_, 2 * size_phi_, 3 * size_phi_)
    , u_(u_space_, u_view_, "u")
    , p_(p_space_, p_view_, "p")
    , Re_(Re)
    , Fa_inv_(1. / Fa)
    , xi_(xi)
    , kappa_(kappa)
    , c_1_(c_1)
    , Pa_(Pa)
    , beta_(beta)
    , gamma_(gamma)
    , Be_(Be)
    , Ca_(Ca)
    , epsilon_(epsilon)
    , In_(In)
    , vol_domain_((upper_right_[0] - lower_left_[0]) * (upper_right_[1] - lower_left_[1]))
    , S_stokes_(size_u_ + size_p_, size_u_ + size_p_, create_stokes_pattern(u_space_, p_space_))
    , A_stokes_(S_stokes_, 0, size_u_, 0, size_u_)
    , B_stokes_(S_stokes_, 0, size_u_, size_u_, size_u_ + size_p_)
    , BT_stokes_(S_stokes_, size_u_, size_u_ + size_p_, 0, size_u_)
    , M_p_stokes_(size_p_, size_p_, make_element_sparsity_pattern(p_space_, p_space_, grid_view_))
    , stokes_rhs_vector_(size_u_ + size_p_, 0.)
    , stokes_f_vector_(stokes_rhs_vector_, 0, size_u_)
    , stokes_g_vector_(stokes_rhs_vector_, size_u_, size_u_ + size_p_)
    , p_basis_integrated_vector_(size_p_)
    , u_dirichlet_constraints_(make_dirichlet_constraints(u_space_, boundary_info_))
    , phi_dirichlet_constraints_(make_dirichlet_constraints(phi_space_, boundary_info_))
    , A_stokes_operator_(
          std::make_shared<MatrixOperator<MatrixViewType, PGV, d>>(grid_view_, u_space_, u_space_, A_stokes_))
    , ofield_submatrix_pattern_(make_element_sparsity_pattern(u_space_, u_space_, grid_view_))
    , S_ofield_(2 * size_u_, 2 * size_u_, create_ofield_pattern(size_u_, ofield_submatrix_pattern_), 100)
    , M_ofield_(size_u_, size_u_, ofield_submatrix_pattern_)
    , C_ofield_elliptic_part_(size_u_, size_u_, ofield_submatrix_pattern_)
    , C_ofield_linear_part_(size_u_, size_u_, ofield_submatrix_pattern_)
    , S_ofield_00_(S_ofield_, 0, size_u_, 0, size_u_)
    , S_ofield_01_(S_ofield_, 0, size_u_, size_u_, 2 * size_u_)
    , S_ofield_10_(S_ofield_, size_u_, 2 * size_u_, 0, size_u_)
    , S_ofield_11_(S_ofield_, size_u_, 2 * size_u_, size_u_, 2 * size_u_)
    , S_ofield_00_operator_(
          std::make_shared<MatrixOperator<MatrixViewType, PGV, d>>(grid_view_, u_space_, u_space_, S_ofield_00_))
    , S_ofield_10_operator_(
          std::make_shared<MatrixOperator<MatrixViewType, PGV, d>>(grid_view_, u_space_, u_space_, S_ofield_10_))
    , M_ofield_operator_(
          std::make_shared<MatrixOperator<MatrixType, PGV, d>>(grid_view_, u_space_, u_space_, M_ofield_))
    , ofield_rhs_vector_(2 * size_u_, 0., 100)
    , ofield_old_result_(num_cells_, EigenVectorType(2 * size_u_, 0.))
    , ofield_f_vector_(ofield_rhs_vector_, 0, size_u_)
    , ofield_g_vector_(ofield_rhs_vector_, size_u_, 2 * size_u_)
    , stokes_solver_(std::make_shared<StokesSolverType>())
    , ofield_solver_(std::make_shared<SolverType>())
    , pfield_solver_(std::make_shared<SolverType>())
    , linearize_(linearize)
    , pfield_submatrix_pattern_(make_element_sparsity_pattern(phi_space_, phi_space_, grid_view_))
    , S_pfield_(3 * size_phi_, 3 * size_phi_, create_pfield_pattern(size_phi_, pfield_submatrix_pattern_), 100)
    , M_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_)
    , E_pfield_(size_phi_, size_phi_, pfield_submatrix_pattern_)
    , A_pfield_linear_part_(size_phi_, size_phi_, pfield_submatrix_pattern_)
    , A_pfield_nonlinear_part_(size_phi_, size_phi_, pfield_submatrix_pattern_)
    , J_pfield_linear_part_(size_phi_, size_phi_, pfield_submatrix_pattern_)
    , S_pfield_00_(S_pfield_, 0, size_phi_, 0, size_phi_)
    , S_pfield_01_(S_pfield_, 0, size_phi_, size_phi_, 2 * size_phi_)
    , S_pfield_10_(S_pfield_, size_phi_, 2 * size_phi_, 0, size_phi_)
    , S_pfield_11_(S_pfield_, size_phi_, 2 * size_phi_, size_phi_, 2 * size_phi_)
    , S_pfield_12_(S_pfield_, size_phi_, 2 * size_phi_, 2 * size_phi_, 3 * size_phi_)
    , S_pfield_20_(S_pfield_, 2 * size_phi_, 3 * size_phi_, 0, size_phi_)
    , S_pfield_22_(S_pfield_, 2 * size_phi_, 3 * size_phi_, 2 * size_phi_, 3 * size_phi_)
    , S_pfield_00_operator_(
          std::make_shared<MatrixOperator<MatrixViewType, PGV, 1>>(grid_view_, phi_space_, phi_space_, S_pfield_00_))
    , S_pfield_10_operator_(
          std::make_shared<MatrixOperator<MatrixViewType, PGV, 1>>(grid_view_, phi_space_, phi_space_, S_pfield_10_))
    , A_pfield_nonlinear_part_operator_(std::make_shared<MatrixOperator<MatrixType, PGV, 1>>(
          grid_view_, phi_space_, phi_space_, A_pfield_nonlinear_part_))
    , pfield_rhs_vector_(3 * size_phi_, 0., 100)
    , pfield_old_result_(num_cells_, EigenVectorType(3 * size_phi_, 0.))
    , pfield_f_vector_(pfield_rhs_vector_, 2 * size_phi_, 3 * size_phi_)
    , pfield_g_vector_(pfield_rhs_vector_, 0, size_phi_)
    , pfield_h_vector_(pfield_rhs_vector_, size_phi_, 2 * size_phi_)
    , u_discr_func_(u_space_)
    , u_discr_func_local_(std::make_shared<StokesPerThread>())
    , P_discr_func_local_(std::make_shared<OfieldPerThread>(num_cells_))
    , Pnat_discr_func_local_(std::make_shared<OfieldPerThread>(num_cells_))
    , phi_discr_func_local_(std::make_shared<PfieldPerThread>(num_cells_))
    , phinat_discr_func_local_(std::make_shared<PfieldPerThread>(num_cells_))
    , mu_discr_func_local_(std::make_shared<PfieldPerThread>(num_cells_))
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

    /*************************************************************************************************
     ******************************* create variables, set initial values ****************************
     *************************************************************************************************/

    // On the non-periodic boundaries, use Dirichlet boundary conditions u = 0 and \phi = -1, Neumann boundary
    // conditions for the other variables
    const XT::Functions::ConstantFunction<d> minus_one(-1.);
    XT::Grid::AllDirichletBoundaryInfo<PI> all_dirichlet_boundary_info;
    default_interpolation(*u_initial_func, u_);
    for (size_t kk = 0; kk < num_cells_; kk++) {
      P_.emplace_back(make_discrete_function(u_space_, P_view_, "P_" + XT::Common::to_string(kk)));
      Pnat_.emplace_back(make_discrete_function(u_space_, Pnat_view_, "Pnat_" + XT::Common::to_string(kk)));
      phi_.emplace_back(make_discrete_function(phi_space_, phi_view_, "phi_" + XT::Common::to_string(kk)));
      phinat_.emplace_back(make_discrete_function(phi_space_, phinat_view_, "phinat_" + XT::Common::to_string(kk)));
      mu_.emplace_back(make_discrete_function(phi_space_, mu_view_, "mu_" + XT::Common::to_string(kk)));
      P_discr_func_.emplace_back(u_space_);
      Pnat_discr_func_.emplace_back(u_space_);
      phi_discr_func_.emplace_back(phi_space_);
      phinat_discr_func_.emplace_back(phi_space_);
      mu_discr_func_.emplace_back(phi_space_);
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

    MatrixOperator<MatrixViewType, PGV, 1, 1, d> B_stokes_operator(grid_view_, p_space_, u_space_, B_stokes_);
    MatrixOperator<MatrixType, PGV, 1, 1, 1> M_p_stokes_operator(grid_view_, p_space_, p_space_, M_p_stokes_);
    // calculate A_{ij} as \int \nabla v_i \nabla v_j
    // A_stokes_operator_->append(LocalElementIntegralBilinearForm<E, d>(LocalSymmetricEllipticIntegrand<E>(1.)));
    A_stokes_operator_->append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>()));
    // calculate B_{ij} as \int \nabla p_i div(v_j)
    B_stokes_operator.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1>(
        LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
    M_p_stokes_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>()));

    auto p_basis_integrated_functional = make_vector_functional(p_space_, p_basis_integrated_vector_);
    const XT::Functions::ConstantGridFunction<E> one_function(1);
    p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), one_function)));
    B_stokes_operator.append(p_basis_integrated_functional);

    // Dirichlet constrainst for u
    A_stokes_operator_->append(u_dirichlet_constraints_);
    // assemble everything
    A_stokes_operator_->assemble(true);
    B_stokes_operator.assemble(true);
    M_p_stokes_operator.assemble(true);

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
    M_ofield_operator_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(1.)));
    M_ofield_operator_->assemble(true);
    // set S_11 = D = M
    S_ofield_11_ = M_ofield_;
    // calculate S_01 = B
    S_ofield_01_ = M_ofield_;
    S_ofield_01_ *= 1. / kappa_;

    MatrixOperator<MatrixType, PGV, d> ofield_elliptic_operator(
        grid_view_, u_space_, u_space_, C_ofield_elliptic_part_);
    ofield_elliptic_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(-1. / Pa_)));
    ofield_elliptic_operator.assemble(true);
    ofield_solver_->analyzePattern(S_ofield_.backend());

    /*************************************************************************************************
     **************************************** Phasefield *********************************************
     *************************************************************************************************/

    MatrixOperator<MatrixType, PGV, 1> M_pfield_operator(grid_view_, phi_space_, phi_space_, M_pfield_);
    M_pfield_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(1.)));
    MatrixOperator<MatrixType, PGV, 1> E_pfield_operator(grid_view_, phi_space_, phi_space_, E_pfield_);
    E_pfield_operator.append(LocalElementIntegralBilinearForm<E, 1>(LocalLaplaceIntegrand<E, 1>(1.)));
    M_pfield_operator.append(phi_dirichlet_constraints_);
    M_pfield_operator.assemble(true);
    E_pfield_operator.assemble(true);
    A_pfield_linear_part_ = E_pfield_ * epsilon_;
    J_pfield_linear_part_ = E_pfield_ * 1. / Be_;
    J_pfield_linear_part_ += M_pfield_ * 1. / Ca_;
    // Set matrix S_{22} = C = M
    S_pfield_22_ = M_pfield_;
    // Set matrix S_{11} = H = M
    S_pfield_11_ = M_pfield_;
    // Set matrix S_{01} = E
    S_pfield_01_ = E_pfield_;
    S_pfield_01_ *= gamma_;

    // apply Dirichlet constraints to linear part
    for (const auto& DoF : phi_dirichlet_constraints_.dirichlet_DoFs()) {
      A_pfield_linear_part_.clear_row(DoF);
      S_pfield_22_.clear_row(DoF);
    }

    pfield_solver_->analyzePattern(S_pfield_.backend());
  } // constructor

  static XT::Common::FieldVector<R, d> get_lower_left(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return {{0., 0.}};
    else if (testcase == "two_cells")
      return {{0., 0.}};
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return FieldVector<R, d>();
  }

  static XT::Common::FieldVector<R, d> get_upper_right(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return {{160., 40.}};
    else if (testcase == "two_cells")
      return {{50., 50.}};
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return FieldVector<R, d>();
  }

  static std::string get_periodic_directions(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return "01";
    else if (testcase == "two_cells")
      return "00";
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return "";
  }

  static size_t get_num_cells(const std::string& testcase)
  {
    if (testcase == "single_cell")
      return 1;
    else if (testcase == "two_cells")
      return 2;
    else
      DUNE_THROW(Dune::NotImplemented, "Unknown testcase");
    return 0;
  }

  size_t num_cells() const
  {
    return num_cells_;
  }

  std::vector<std::vector<VectorType>> solve(const double dt,
                                             const bool write,
                                             const double write_step,
                                             const std::string filename = "cellmodel",
                                             const bool subsampling = true)
  {
    std::vector<std::vector<VectorType>> ret(3);
    ret[0].push_back(pfield_vector_);
    ret[1].push_back(ofield_vector_);
    ret[2].push_back(stokes_vector_);
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
        ret[0].push_back(solve_pfield(ret[0].back(), kk));
        set_pfield_variables(kk, ret[0].back());
        std::cout << "Pfield " << kk << " done" << std::endl;
        prepare_ofield_operator(dt, kk);
        ret[1].push_back(solve_ofield(ret[1].back(), kk));
        set_ofield_variables(kk, ret[1].back());
        std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      prepare_stokes_operator();
      ret[2].push_back(solve_stokes());
      set_stokes_variables(ret[2].back());
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

  std::vector<std::vector<VectorType>> next_n_timesteps(const size_t n, const double dt)
  {
    std::vector<std::vector<VectorType>> ret(3);
    size_t count = 0;
    if (XT::Common::is_zero(t_)) {
      ret[0].push_back(pfield_vector_);
      ret[1].push_back(ofield_vector_);
      ret[2].push_back(stokes_vector_);
      // Hack to avoid adding initial_values twice
      ++count;
      t_ = 1e-100;
    }
    // Undo hack to avoid adding initial_values twice
    if (XT::Common::is_zero(t_ - 1e-100))
      t_ = 0.;

    // implicit Euler timestepping
    assert(Dune::XT::Common::FloatCmp::ge(t_end_, t_));

    while (Dune::XT::Common::FloatCmp::lt(t_, t_end_) && count < n) {
      double max_dt = dt;
      // match saving times and t_end_ exactly
      if (Dune::XT::Common::FloatCmp::gt(t_ + dt, t_end_))
        max_dt = t_end_ - t_;
      double actual_dt = std::min(dt, max_dt);

      // do a timestep
      for (size_t kk = 0; kk < num_cells_; ++kk) {
        prepare_pfield_operator(dt, kk);
        ret[0].push_back(solve_pfield(ret[0].back(), kk));
        set_pfield_variables(kk, ret[0].back());
        std::cout << "Pfield " << kk << " done" << std::endl;
        prepare_ofield_operator(dt, kk);
        ret[1].push_back(solve_ofield(ret[1].back(), kk));
        set_ofield_variables(kk, ret[1].back());
        std::cout << "Ofield " << kk << " done" << std::endl;
      }

      // stokes system
      prepare_stokes_operator();
      ret[2].push_back(solve_stokes());
      set_stokes_variables(ret[2].back());
      std::cout << "Stokes done" << std::endl;

      ++count;
      t_ += actual_dt;
    } // while (t_ < t_end_)
    pfield_vector_ = ret[0].back();
    ofield_vector_ = ret[1].back();
    stokes_vector_ = ret[2].back();
    return ret;
  }

  // applies the pfield mass matrix to phi, phinat, mu
  // To calculate the sum of the squared L2 products of phi, phinat and mu, calculate the inner product of the result
  // with vec.
  VectorType apply_pfield_product_operator(const VectorType& vec) const
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
  VectorType apply_ofield_product_operator(const VectorType& vec) const
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
  VectorType apply_stokes_product_operator(const VectorType& vec) const
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

  bool linear() const
  {
    return linearize_;
  }

  bool finished() const
  {
    return XT::Common::FloatCmp::eq(t_end_, t_);
  }

  void visualize(const std::string& prefix,
                 const size_t step,
                 const double t,
                 const bool subsampling = true,
                 const bool vtu = true,
                 const bool txt = false) const
  {
    auto vtk_writer = u_.create_vtkwriter(u_.space().grid_view(), subsampling);
    std::string postfix = "_" + XT::Common::to_string(step);
    if (vtu) {
      u_.add_to_vtkwriter(*vtk_writer);
      p_.add_to_vtkwriter(*vtk_writer);
      for (size_t kk = 0; kk < num_cells_; ++kk) {
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

  VectorType stokes_vector()
  {
    return stokes_vector_;
  }

  VectorType ofield_vector()
  {
    return ofield_vector_;
  }

  VectorType pfield_vector()
  {
    return pfield_vector_;
  }

  static XT::LA::SparsityPatternDefault create_ofield_pattern(const size_t n,
                                                              const XT::LA::SparsityPatternDefault& submatrix_pattern)
  {
    XT::LA::SparsityPatternDefault pattern(2 * n);
    for (size_t ii = 0; ii < n; ++ii)
      for (const auto& jj : submatrix_pattern.inner(ii)) {
        pattern.insert(ii, jj);
        pattern.insert(ii, n + jj);
        pattern.insert(n + ii, jj);
        pattern.insert(n + ii, n + jj);
      }
    pattern.sort();
    return pattern;
  }

  static XT::LA::SparsityPatternDefault create_pfield_pattern(const size_t n,
                                                              const XT::LA::SparsityPatternDefault& submatrix_pattern)
  {
    // Use same pattern for all submatrices
    XT::LA::SparsityPatternDefault pattern(3 * n);
    for (size_t ii = 0; ii < n; ++ii)
      for (const auto& jj : submatrix_pattern.inner(ii)) {
        pattern.insert(ii, jj); // S_{00}
        pattern.insert(ii, n + jj); // S_{01}
        pattern.insert(n + ii, jj); // S_{10}
        pattern.insert(n + ii, n + jj); // S_{11}
        pattern.insert(n + ii, 2 * n + jj); // S_{12}
        pattern.insert(2 * n + ii, jj); // S_{20}
        pattern.insert(2 * n + ii, 2 * n + jj); // S_{22}
      }
    pattern.sort();
    return pattern;
  }

  static XT::LA::SparsityPatternDefault create_stokes_pattern(const SpaceInterface<PGV, d, 1, R>& u_space,
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

  void set_pfield_variables(const size_t ll, const VectorType& pfield_vec)
  {
    DUNE_THROW_IF(ll >= num_cells_, InvalidStateException, "");
    for (size_t nn = 0; nn < size_phi_; nn++) {
      phi_[ll].dofs().vector().set_entry(nn, pfield_vec.get_entry(nn));
      phinat_[ll].dofs().vector().set_entry(nn, pfield_vec.get_entry(size_phi_ + nn));
      mu_[ll].dofs().vector().set_entry(nn, pfield_vec.get_entry(2 * size_phi_ + nn));
    }
  }

  void set_ofield_variables(const size_t ll, const VectorType& ofield_vec)
  {
    DUNE_THROW_IF(ll >= num_cells_, InvalidStateException, "");
    for (size_t nn = 0; nn < size_u_; nn++) {
      P_[ll].dofs().vector().set_entry(nn, ofield_vec.get_entry(nn));
      Pnat_[ll].dofs().vector().set_entry(nn, ofield_vec.get_entry(size_u_ + nn));
    }
  }

  void set_stokes_variables(const VectorType& stokes_vec)
  {
    for (size_t ii = 0; ii < size_u_; ++ii)
      u_.dofs().vector()[ii] = stokes_vec[ii];
    for (size_t ii = 0; ii < size_p_; ++ii)
      p_.dofs().vector()[ii] = stokes_vec[size_u_ + ii];
  }

  void prepare_stokes_operator()
  {
    auto begin = std::chrono::steady_clock::now();
    auto f_functional = make_vector_functional(u_space_, stokes_f_vector_);

    u_discr_func_.dofs().vector() = u_.dofs().vector();
    for (size_t kk = 0; kk < num_cells_; kk++) {
      phi_discr_func_[kk].dofs().vector() = phi_[kk].dofs().vector();
      phinat_discr_func_[kk].dofs().vector() = phinat_[kk].dofs().vector();
      P_discr_func_[kk].dofs().vector() = P_[kk].dofs().vector();
      Pnat_discr_func_[kk].dofs().vector() = Pnat_[kk].dofs().vector();
    }

    // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
    f_functional.append(LocalElementIntegralFunctional<E, d>(
        /*order*/ [& u_space_ = P_[0].space()](
                      const auto& test_basis,
                      const auto& param) { return 3 * u_space_.max_polorder() + test_basis.order(param); },
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
    A_stokes_operator_->clear();
    A_stokes_operator_->append(f_functional);
    stokes_f_vector_ *= 0.;
    A_stokes_operator_->assemble(true);
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Assembling Stokes took: " << time.count() << " s!" << std::endl;

    // apply dirichlet constraints for u.
    u_dirichlet_constraints_.apply(stokes_f_vector_);
  }

  VectorType apply_stokes_operator(VectorType source) const
  {
    VectorType ret(size_u_ + size_p_, 0.);
    S_stokes_.mv(source, ret);
    ret -= stokes_rhs_vector_;
    return ret;
  }

  VectorType solve_stokes() const
  {
    // now solve the system
    auto begin = std::chrono::steady_clock::now();
    EigenVectorType ret(size_u_ + size_p_);
    ret.backend() = stokes_solver_->solve(stokes_rhs_vector_.backend());
    std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
    std::cout << "Solving Stokes took: " << time.count() << " s!" << std::endl;

    // ensure int_\Omega p = 0 (TODO: remove, not necessary as p is not used anywhere)
    // auto p_integral = p_basis_integrated_vector_ * p_.dofs().vector();
    // auto p_correction = make_discrete_function<VectorType>(p_space_, "p_corr");
    // XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain_);
    // default_interpolation(const_p_integral_func, p_correction);
    // p_ -= p_correction;

    return XT::Common::convert_to<VectorType>(ret);
  }

  void prepare_ofield_operator(const double dt, const size_t ll)
  {
    u_discr_func_.dofs().vector() = u_.dofs().vector();
    P_discr_func_[ll].dofs().vector() = P_[ll].dofs().vector();
    phi_discr_func_[ll].dofs().vector() = phi_[ll].dofs().vector();
    assemble_ofield_rhs(dt, ll);
    assemble_ofield_linear_jacobian(dt, ll);
  }

  VectorType apply_ofield_operator(const VectorType& source, const size_t ll) const
  {
    // linear part
    VectorType residual(source.size());
    S_ofield_.mv(source, residual);
    // subtract rhs
    residual -= ofield_rhs_vector_;
    if (!linearize_) {
      fill_tmp_ofield(ll, source);
      // nonlinear part
      VectorViewType res0_vec(residual, 0, size_u_);
      VectorViewType res1_vec(residual, size_u_, 2 * size_u_);
      const auto res0 = make_discrete_function(u_space_, res0_vec);
      const auto res1 = make_discrete_function(u_space_, res1_vec);
      auto nonlinear_res_functional = make_vector_functional(u_space_, res1_vec);
      XT::Functions::GenericGridFunction<E, d, 1> nonlinear_res_pf(
          /*order = */ 3 * u_space_.max_polorder(),
          /*post_bind_func*/
          [ll, this](const E& element) { this->bind_P(ll, element); },
          /*evaluate_func*/
          [ll, factor = -c_1_ * 1. / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
            // evaluate P, divP
            const auto P_n = this->eval_P(ll, x_local, param);
            return P_n * (factor * (P_n * P_n));
          });
      nonlinear_res_functional.append(LocalElementIntegralFunctional<E, d>(
          local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), nonlinear_res_pf)));
      S_ofield_00_operator_->clear();
      S_ofield_00_operator_->append(nonlinear_res_functional);
      S_ofield_00_operator_->assemble(true);
    }
    // relative error if l2_norm is > 1, else absolute error
    return residual;
  }

  void assemble_ofield_rhs(const double dt, const size_t ll)
  {
    M_ofield_.mv(P_[ll].dofs().vector(), ofield_f_vector_);
    ofield_f_vector_ /= dt;

    auto g_functional = make_vector_functional(u_space_, ofield_g_vector_);
    ofield_g_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, d> g(
        /*order = */ 3 * u_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) {
          this->bind_phi(ll, element);
          if (this->linearize_) {
            this->bind_P(ll, element);
          }
        },
        /*evaluate_func*/
        [ll, factor1 = beta_ * 1. / Pa_, factor2 = -2. * c_1_ / Pa_, this](const DomainType& x_local,
                                                                           const XT::Common::Parameter& param) {
          // evaluate rhs terms
          const auto grad_phi = this->grad_phi(ll, x_local, param);
          auto ret = grad_phi;
          ret *= factor1;
          if (linearize_) {
            const auto P_n = this->eval_P(ll, x_local, param);
            auto ret2 = P_n;
            ret2 *= factor2 * (P_n * P_n);
            ret += ret2;
          }
          return ret;
        });
    g_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), g)));
    S_ofield_10_operator_->clear();
    S_ofield_10_operator_->append(g_functional);
    S_ofield_10_operator_->assemble(true);
  }

  void assemble_ofield_linear_jacobian(const double dt, const size_t ll)
  {
    // assemble matrix S_{00} = M/dt + A
    S_ofield_00_operator_->clear();
    S_ofield_00_ = M_ofield_;
    S_ofield_00_ *= 1 / dt;
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
    S_ofield_00_operator_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(Omega_minus_xi_D_transposed)));
    S_ofield_00_operator_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementGradientValueIntegrand<E, d>(u_)));
    S_ofield_00_operator_->assemble(true);

    // calculate linear part S_10 = C
    S_ofield_10_operator_->clear();
    S_ofield_10_ = C_ofield_elliptic_part_;
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_inv_phi(
        /*order = */ u_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) { this->bind_phi(ll, element); },
        /*evaluate_func*/
        [ll, factor = c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto phi = this->eval_phi(ll, x_local, param);
          return factor * phi;
        });
    S_ofield_10_operator_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(c1_Pa_inv_phi)));
    S_ofield_10_operator_->assemble(true);
    C_ofield_linear_part_ = S_ofield_10_;

    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_ofield_nonlinear_jacobian(ofield_vector(), ll);
  }

  void assemble_ofield_nonlinear_jacobian(const VectorType& source, const size_t ll) const
  {
    S_ofield_10_operator_->clear();
    ConstVectorViewType P_vec(source, 0, size_u_);
    const auto P = make_discrete_function(u_space_, P_vec);
    XT::Functions::GenericGridFunction<E, 1, 1> c1_Pa_P2(
        /*order = */ 2. * u_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) { this->bind_P(ll, element); },
        /*evaluate_func*/
        [ll, factor = -c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = this->eval_P(ll, x_local, param);
          return factor * P_n.two_norm2();
        });
    S_ofield_10_operator_->append(LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(c1_Pa_P2)));
    XT::Functions::GenericGridFunction<E, d, d> minus_two_frac_c1_Pa_Pn_otimes_Pn(
        /*order = */ 2 * u_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) { this->bind_P(ll, element); },
        /*evaluate_func*/
        [ll, factor = -2. * c_1_ / Pa_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const auto P_n = this->eval_P(ll, x_local, param);
          FieldMatrix<R, d, d> ret;
          for (size_t ii = 0; ii < d; ++ii)
            for (size_t jj = 0; jj < d; ++jj)
              ret[ii][jj] = P_n[ii] * P_n[jj];
          ret *= factor;
          return ret;
        });
    // Pn_otimes_Pn is symmetric, so no need to transpose
    S_ofield_10_operator_->append(
        LocalElementIntegralBilinearForm<E, d>(LocalElementProductIntegrand<E, d>(minus_two_frac_c1_Pa_Pn_otimes_Pn)));
    S_ofield_10_operator_->assemble(true);
  }

  void revert_ofield_jacobian_to_linear() const
  {
    S_ofield_10_ = C_ofield_linear_part_;
  }

  VectorType solve_ofield_linear_system(const VectorType& rhs, const size_t ll) const
  {
    //    std::ofstream S_file("S_" + XT::Common::to_string(dt) + ".txt");
    //    S_file << S_ << std::endl;
    //    S_file.close();
    //    DUNE_THROW(NotImplemented, "");
    //    const auto ret = XT::LA::solve(S_, rhs_vector_, XT::LA::SolverOptions<MatrixType>::options("lu.umfpack"));
    //      ofield_update_ = XT::LA::solve(S_, ofield_residual_);
    EigenVectorType update(rhs.size());
    ofield_solver_->compute(S_ofield_.backend());
    const auto rhs_eig = XT::Common::convert_to<EigenVectorType>(rhs);
    update.backend() = ofield_solver_->solveWithGuess(rhs_eig.backend(), ofield_old_result_[ll].backend());
    ofield_old_result_[ll] = update;
    return XT::Common::convert_to<VectorType>(update);
  }

  double ofield_residual_norm(const VectorType& residual, double l2_ref_P, double l2_ref_Pnat) const
  {
    l2_ref_P = l2_ref_P < 1. ? 1. : l2_ref_P;
    l2_ref_Pnat = l2_ref_Pnat < 1. ? 1. : l2_ref_Pnat;
    ConstVectorViewType res0_vec(residual, 0, size_u_);
    ConstVectorViewType res1_vec(residual, size_u_, 2 * size_u_);
    const auto res0 = make_discrete_function(u_space_, res0_vec);
    const auto res1 = make_discrete_function(u_space_, res1_vec);
    return l2_norm(grid_view_, res0) / l2_ref_P + l2_norm(grid_view_, res1) / l2_ref_Pnat;
  }

  VectorType solve_ofield(const VectorType& source, const size_t ll)
  {
    if (linearize_) {
      auto residual = apply_ofield_operator(source, ll);
      residual *= -1;
      const auto update = solve_ofield_linear_system(residual, ll);
      return source + update;
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 1000;
      const auto max_dampening_iter = 1000;

      auto l2_norm_P = l2_norm(grid_view_, P_[ll]);
      auto l2_norm_Pnat = l2_norm(grid_view_, Pnat_[ll]);

      // ********* compute residual *********
      auto begin = std::chrono::steady_clock::now();
      auto residual = apply_ofield_operator(source, ll);
      auto res_norm = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
      std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
      std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

      size_t iter = 0;
      VectorType x_n = source;
      VectorType x_n_plus_1 = source;
      VectorType update;
      while (true) {
        if (res_norm < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_ofield_nonlinear_jacobian(x_n, ll);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        residual *= -1.;
        update = solve_ofield_linear_system(residual, ll);

        DUNE_THROW_IF(
            iter >= max_iter, Exceptions::operator_error, "max iterations reached!\n|residual|_l2 = " << res_norm);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        revert_ofield_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res_norm) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res_norm << "\nl = " << iter << "\n");
          x_n_plus_1 = x_n + update * lambda;
          residual = apply_ofield_operator(x_n_plus_1, ll);
          candidate_res = ofield_residual_norm(residual, l2_norm_P, l2_norm_Pnat);
          std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        x_n = x_n_plus_1;
        res_norm = candidate_res;
        std::cout << "Current res: " << candidate_res << std::endl;
        iter += 1;
      } // while (true)
      return x_n;
    }
  }

  double
  pfield_residual_norm(const VectorType& residual, double l2_ref_phi, double l2_ref_phinat, double l2_ref_mu) const
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


  void prepare_pfield_operator(const double dt, const size_t ll)
  {
    u_discr_func_.dofs().vector() = u_.dofs().vector();
    for (size_t kk = 0; kk < num_cells_; kk++) {
      phi_discr_func_[kk].dofs().vector() = phi_[kk].dofs().vector();
      mu_discr_func_[kk].dofs().vector() = mu_[kk].dofs().vector();
    }
    assemble_pfield_rhs(dt, ll);
    assemble_pfield_linear_jacobian(dt, ll);
  }

  VectorType apply_pfield_operator(const VectorType& source, const size_t ll)
  {
    // linear part
    VectorType residual(source.size());
    S_pfield_.mv(source, residual);
    // subtract rhs
    residual -= pfield_rhs_vector_;
    VectorViewType res2_vec(residual, 2 * size_phi_, 3 * size_phi_);

    if (!linearize_) {
      fill_tmp_pfield(ll, source);
      // nonlinear part
      VectorViewType res0_vec(residual, 0, size_phi_);
      VectorViewType res1_vec(residual, size_phi_, 2 * size_phi_);
      const auto res0 = make_discrete_function(phi_space_, res0_vec);
      const auto res1 = make_discrete_function(phi_space_, res1_vec);
      const auto res2 = make_discrete_function(phi_space_, res2_vec);
      auto nonlinear_res1_functional = make_vector_functional(phi_space_, res1_vec);
      auto nonlinear_res2_functional = make_vector_functional(phi_space_, res2_vec);
      const auto Bfunc = [epsilon_inv = 1. / epsilon_,
                          this](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi_n = this->eval_phi(kk, x_local, param);
        return epsilon_inv * std::pow(std::pow(phi_n, 2) - 1, 2);
      };
      const auto wfunc = [this](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
        const auto phi_n = this->eval_phi(kk, x_local, param);
        if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
          return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
        else
          return 0.;
      };
      XT::Functions::GenericGridFunction<E, 1, 1> nonlinear_res_pf1(
          /*order = */ 3 * phi_space_.max_polorder(),
          /*post_bind_func*/
          [ll, this](const E& element) {
            this->bind_phi(ll, element);
            this->bind_mu(ll, element);
            if (this->num_cells_ > 1) {
              for (size_t kk = 0; kk < this->num_cells_; ++kk)
                this->bind_phi(kk, element);
            }
          },
          /*evaluate_func*/
          [wfunc,
           Bfunc,
           ll,
           In_inv = 1. / In_,
           eps_inv = 1. / epsilon_,
           num_cells = num_cells_,
           inv_Be_eps2 = 1. / (Be_ * std::pow(epsilon_, 2)),
           this](const DomainType& x_local, const XT::Common::Parameter& param) {
            // evaluate P, divP
            const auto phi_n = this->eval_phi(ll, x_local, param);
            const auto mu_n = this->eval_mu(ll, x_local, param);
            auto ret = inv_Be_eps2 * (3. * phi_n * phi_n - 1) * mu_n;
            if (num_cells > 1) {
              R wsum = 0.;
              R Bsum = 0.;
              for (size_t kk = 0; kk < num_cells; ++kk) {
                if (kk != ll) {
                  wsum += wfunc(kk, x_local, param);
                  Bsum += Bfunc(kk, x_local, param);
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
          [ll, this](const E& element) { this->bind_phi(ll, element); },
          /*evaluate_func*/
          [ll, inv_eps = 1. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
            // evaluate P, divP
            const auto phi_n = this->eval_phi(ll, x_local, param);
            return inv_eps * (phi_n * phi_n - 1) * phi_n;
          });
      nonlinear_res1_functional.append(LocalElementIntegralFunctional<E, 1>(
          local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), nonlinear_res_pf1)));
      nonlinear_res2_functional.append(LocalElementIntegralFunctional<E, 1>(
          local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), nonlinear_res_pf2)));
      S_pfield_00_operator_->clear();
      S_pfield_00_operator_->append(nonlinear_res1_functional);
      S_pfield_00_operator_->append(nonlinear_res2_functional);
      S_pfield_00_operator_->assemble(true);
    }
    phi_dirichlet_constraints_.apply(res2_vec);
    // relative error if l2_norm is > 1, else absolute error
    return residual;
  }

  void fill_tmp_ofield(const size_t ll, const VectorType& source) const
  {
    ConstVectorViewType P_vec(source, 0, size_u_);
    P_discr_func_[ll].dofs().vector() = P_vec;
  }

  void fill_tmp_pfield(const size_t ll, const VectorType& source) const
  {
    ConstVectorViewType phi_vec(source, 0, size_phi_);
    ConstVectorViewType mu_vec(source, 2 * size_phi_, 3 * size_phi_);
    phi_discr_func_[ll].dofs().vector() = phi_vec;
    mu_discr_func_[ll].dofs().vector() = mu_vec;
  }

  void bind_u(const E& element) const
  {
    auto& u_local = **u_discr_func_local_;
    if (!u_local)
      u_local = u_discr_func_.local_function();
    u_local->bind(element);
  }

  DomainRetType eval_u(const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    return (**u_discr_func_local_)->evaluate(x_local, param);
  }

  JacobianRetType grad_u(const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    return (**u_discr_func_local_)->jacobian(x_local, param);
  }

  void bind_P(const size_t ll, const E& element) const
  {
    auto& P_local_ll = (**P_discr_func_local_)[ll];
    if (!P_local_ll)
      P_local_ll = P_discr_func_[ll].local_function();
    P_local_ll->bind(element);
  }

  DomainRetType eval_P(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& P_local_ll = (**P_discr_func_local_)[ll];
    return P_local_ll->evaluate(x_local, param);
  }

  JacobianRetType grad_P(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& P_local_ll = (**P_discr_func_local_)[ll];
    return P_local_ll->jacobian(x_local, param);
  }

  void bind_Pnat(const size_t ll, const E& element) const
  {
    auto& Pnat_local_ll = (**Pnat_discr_func_local_)[ll];
    if (!Pnat_local_ll)
      Pnat_local_ll = Pnat_discr_func_[ll].local_function();
    Pnat_local_ll->bind(element);
  }

  DomainRetType eval_Pnat(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param) const
  {
    auto& Pnat_local_ll = (**Pnat_discr_func_local_)[ll];
    return Pnat_local_ll->evaluate(x_local, param);
  }

  void bind_phi(const size_t ll, const E& element) const
  {
    auto& phi_local_ll = (**phi_discr_func_local_)[ll];
    if (!phi_local_ll)
      phi_local_ll = phi_discr_func_[ll].local_function();
    phi_local_ll->bind(element);
  }

  R eval_phi(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phi_local_ll = (**phi_discr_func_local_)[ll];
    return phi_local_ll->evaluate(x_local, param)[0];
  }

  DomainRetType grad_phi(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phi_local_ll = (**phi_discr_func_local_)[ll];
    return phi_local_ll->jacobian(x_local, param)[0];
  }

  void bind_phinat(const size_t ll, const E& element) const
  {
    auto& phinat_local_ll = (**phinat_discr_func_local_)[ll];
    if (!phinat_local_ll)
      phinat_local_ll = phinat_discr_func_[ll].local_function();
    phinat_local_ll->bind(element);
  }

  R eval_phinat(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& phinat_local_ll = (**phinat_discr_func_local_)[ll];
    return phinat_local_ll->evaluate(x_local, param)[0];
  }

  void bind_mu(const size_t ll, const E& element) const
  {
    auto& mu_local_ll = (**mu_discr_func_local_)[ll];
    if (!mu_local_ll)
      mu_local_ll = mu_discr_func_[ll].local_function();
    mu_local_ll->bind(element);
  }

  R eval_mu(const size_t ll, const DomainType& x_local, const XT::Common::Parameter& param)
  {
    auto& mu_local_ll = (**mu_discr_func_local_)[ll];
    return mu_local_ll->evaluate(x_local, param)[0];
  }

  void assemble_pfield_rhs(const double dt, const size_t ll)
  {
    S_pfield_00_operator_->clear();
    auto f_functional = make_vector_functional(phi_space_, pfield_f_vector_);
    auto h_functional = make_vector_functional(phi_space_, pfield_h_vector_);

    // calculate f
    if (linearize_)
      pfield_f_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, 1, 1> f_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) { this->bind_phi(ll, element); },
        /*evaluate_func*/
        [ll, two_epsilon_inv = 2. / epsilon_, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate phi_
          const R phi_n = this->eval_phi(ll, x_local, param);
          return two_epsilon_inv * std::pow(phi_n, 3);
        });
    f_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(), f_pf)));
    if (linearize_)
      S_pfield_00_operator_->append(f_functional);

    // calculate g
    M_pfield_.mv(phi_[ll].dofs().vector(), pfield_g_vector_);
    pfield_g_vector_ /= dt;

    // calculate h
    pfield_h_vector_ *= 0.;
    XT::Functions::GenericGridFunction<E, 1, 1> h_pf(
        /*order = */ 3 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) {
          this->bind_P(ll, element);
          if (this->linearize_) {
            this->bind_phi(ll, element);
            this->bind_mu(ll, element);
          }
        },
        /*evaluate_func*/
        [ll, factor0 = 6. / (Be_ * std::pow(epsilon_, 2)), factor1 = -c_1_ / (2. * Pa_), factor2 = -beta_ / Pa_, this](
            const DomainType& x_local, const XT::Common::Parameter& param) {
          // evaluate P, divP
          const auto Pn = this->eval_P(ll, x_local, param);
          const auto grad_P = this->grad_P(ll, x_local, param);
          R div_P(0.);
          for (size_t ii = 0; ii < d; ++ii)
            div_P += grad_P[ii][ii];
          auto ret = factor1 * (Pn * Pn) + factor2 * div_P;
          if (this->linearize_) {
            const auto phi_n = this->eval_phi(ll, x_local, param);
            const auto mu_n = this->eval_mu(ll, x_local, param);
            ret += factor0 * std::pow(phi_n, 2) * mu_n;
          }
          return ret;
        });
    h_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, 1>(1.), h_pf)));
    S_pfield_00_operator_->append(h_functional);

    // assemble rhs
    S_pfield_00_operator_->assemble(true);
  }

  void assemble_pfield_linear_jacobian(const double dt, const size_t ll)
  {
    // assemble matrix S_{00} = M/dt + D
    S_pfield_00_operator_->clear();
    S_pfield_00_ = M_pfield_;
    S_pfield_00_ *= 1. / dt;
    XT::Functions::GenericGridFunction<E, d, 1> minus_u(
        /*order = */ u_space_.max_polorder(),
        /*post_bind_func*/
        [this](const E& element) { this->bind_u(element); },
        /*evaluate_func*/
        [this](const DomainType& x_local, const XT::Common::Parameter& param) {
          auto ret = this->eval_u(x_local, param);
          ret *= -1.;
          return ret;
        });
    S_pfield_00_operator_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementGradientValueIntegrand<E, 1, 1, R, R, R, true>(minus_u)));
    S_pfield_00_operator_->assemble(true);
    // linear part of matrix S_{12} = J
    S_pfield_12_ = J_pfield_linear_part_;
    // linear part of matrix S_{20} = A
    S_pfield_20_ = A_pfield_linear_part_;

    // nonlinear part is equal to linearized part in first iteration
    if (linearize_)
      assemble_pfield_nonlinear_jacobian(pfield_vector(), ll);
  }

  void assemble_pfield_nonlinear_jacobian(const VectorType& source, const size_t ll)
  {
    fill_tmp_pfield(ll, source);
    A_pfield_nonlinear_part_operator_->clear();
    A_pfield_nonlinear_part_ *= 0.;
    XT::Functions::GenericGridFunction<E, 1, 1> A_nonlinear_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) { this->bind_phi(ll, element); },
        /*evaluate_func*/
        [ll, this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const R phi_n = this->eval_phi(ll, x_local, param);
          return (3. * phi_n * phi_n - 1.);
        });
    A_pfield_nonlinear_part_operator_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(A_nonlinear_prefactor)));
    A_pfield_nonlinear_part_operator_->assemble(true);
    A_pfield_nonlinear_part_ *= 1. / epsilon_;
    S_pfield_20_ += A_pfield_nonlinear_part_;
    A_pfield_nonlinear_part_ *= 1. / (Be_ * epsilon_);
    S_pfield_12_ += A_pfield_nonlinear_part_;

    // assemble matrix S_{10} = G
    S_pfield_10_operator_->clear();
    S_pfield_10_ *= 0.;
    const auto Bfunc = [epsilon_inv = 1. / epsilon_,
                        this](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
      const R phi_n = this->eval_phi(kk, x_local, param);
      return epsilon_inv * std::pow(std::pow(phi_n, 2) - 1, 2);
    };
    const auto wfunc = [this](const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param) {
      const R phi_n = this->eval_phi(kk, x_local, param);
      if (XT::Common::FloatCmp::lt(std::abs(phi_n), 1.))
        return std::exp(-0.5 * std::pow(std::log((1 + phi_n) / (1 - phi_n)), 2));
      else
        return 0.;
    };
    XT::Functions::GenericGridFunction<E, 1, 1> G_prefactor(
        /*order = */ 2 * phi_space_.max_polorder(),
        /*post_bind_func*/
        [ll, this](const E& element) {
          this->bind_phi(ll, element);
          this->bind_mu(ll, element);
          if (this->num_cells_ > 1) {
            for (size_t kk = 0; kk < this->num_cells_; ++kk)
              this->bind_phi(kk, element);
          }
        },
        /*evaluate_func*/
        [ll,
         In_inv = 1. / In_,
         eps_inv = 1. / epsilon_,
         six_inv_Be_eps2 = 6. / (Be_ * std::pow(epsilon_, 2)),
         &Bfunc,
         &wfunc,
         this](const DomainType& x_local, const XT::Common::Parameter& param) {
          const R phi_n = this->eval_phi(ll, x_local, param);
          const R mu_n = this->eval_mu(ll, x_local, param);
          auto ret = six_inv_Be_eps2 * phi_n * mu_n;
          if (this->num_cells_ > 1) {
            R wsum = 0.;
            R Bsum = 0.;
            for (size_t kk = 0; kk < this->num_cells_; ++kk) {
              if (kk != ll) {
                wsum += wfunc(kk, x_local, param);
                Bsum += Bfunc(kk, x_local, param);
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
    S_pfield_10_operator_->append(
        LocalElementIntegralBilinearForm<E, 1>(LocalElementProductIntegrand<E, 1>(G_prefactor)));
    S_pfield_10_operator_->assemble(true);

    for (const auto& DoF : phi_dirichlet_constraints_.dirichlet_DoFs())
      S_pfield_20_.unit_row(DoF);
  }

  void revert_pfield_jacobian_to_linear()
  {
    // clear S_{10} = G
    S_pfield_10_ *= 0.;
    // linear part of matrix S_{12} = J
    S_pfield_12_ = J_pfield_linear_part_;
    // linear part of matrix S_{20} = A
    S_pfield_20_ = A_pfield_linear_part_;
  }

  VectorType solve_pfield_linear_system(const VectorType& rhs, const size_t ll)
  {
    EigenVectorType update(rhs.size());
    pfield_solver_->compute(S_pfield_.backend());
    const auto rhs_eig = XT::Common::convert_to<EigenVectorType>(rhs);
    update.backend() = pfield_solver_->solveWithGuess(rhs_eig.backend(), pfield_old_result_[ll].backend());
    pfield_old_result_[ll] = update;
    return XT::Common::convert_to<VectorType>(update);
  }

  VectorType solve_pfield(const VectorType& source, const size_t ll)
  {
    if (linearize_) {
      auto residual = apply_pfield_operator(source, ll);
      residual *= -1;
      const auto update = solve_pfield_linear_system(residual, ll);
      return source + update;
    } else {

      // *********** Newton ******************************
      const auto tol = 1e-10;
      const auto max_iter = 1000;
      const auto max_dampening_iter = 1000;

      const auto l2_norm_phi = l2_norm(grid_view_, phi_[ll]);
      const auto l2_norm_phinat = l2_norm(grid_view_, phinat_[ll]);
      const auto l2_norm_mu = l2_norm(grid_view_, mu_[ll]);

      // ********* compute residual *********
      auto begin = std::chrono::steady_clock::now();
      auto residual = apply_pfield_operator(source, ll);
      auto res_norm = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
      std::chrono::duration<double> time = std::chrono::steady_clock::now() - begin;
      std::cout << "Computing residual took: " << time.count() << " s!" << std::endl;

      size_t iter = 0;
      VectorType x_n = source;
      VectorType x_n_plus_1 = source;
      VectorType update;
      while (true) {
        if (res_norm < tol)
          break;

        // ********** assemble nonlinear part of S = Jacobian ***********
        begin = std::chrono::steady_clock::now();
        assemble_pfield_nonlinear_jacobian(x_n, ll);
        time = std::chrono::steady_clock::now() - begin;
        std::cout << "Assembling nonlinear part of jacobian took: " << time.count() << " s!" << std::endl;

        // *********** solve system *************
        residual *= -1.;
        update = solve_pfield_linear_system(residual, ll);

        DUNE_THROW_IF(
            iter >= max_iter, Exceptions::operator_error, "max iterations reached!\n|residual|_l2 = " << res_norm);

        // apply damping
        size_t k = 0;
        auto candidate_res = 2 * res_norm; // any number such that we enter the while loop at least once
        double lambda = 1;

        // revert jacobian back to linear part to correctly calculate linear part of residual
        revert_pfield_jacobian_to_linear();

        // backtracking line search
        const double gamma = 0.001;
        while (candidate_res > (1 - gamma * lambda) * res_norm) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res_norm << "\nl = " << iter << "\n");
          x_n_plus_1 = x_n + update * lambda;
          residual = apply_pfield_operator(x_n_plus_1, ll);
          candidate_res = pfield_residual_norm(residual, l2_norm_phi, l2_norm_phinat, l2_norm_mu);
          std::cout << "Candidate res: " << candidate_res << std::endl;
          lambda /= 2;
          k += 1;
        }
        x_n = x_n_plus_1;
        res_norm = candidate_res;
        std::cout << "Current res: " << candidate_res << std::endl;
        iter += 1;
      } // while (true)
      return x_n;
    }
  }

  XT::Common::FieldVector<R, d> lower_left_;
  XT::Common::FieldVector<R, d> upper_right_;
  XT::Grid::GridProvider<G> grid_;
  const double t_end_;
  double t_;
  const GV nonperiodic_grid_view_;
  const PGV grid_view_;
  const ContinuousLagrangeSpace<PGV, d, R> u_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> p_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> phi_space_;
  const size_t size_u_;
  const size_t size_p_;
  const size_t size_phi_;
  const size_t num_cells_;
  VectorType stokes_vector_;
  VectorType ofield_vector_;
  VectorType pfield_vector_;
  VectorViewType u_view_;
  VectorViewType p_view_;
  VectorViewType P_view_;
  VectorViewType Pnat_view_;
  VectorViewType phi_view_;
  VectorViewType phinat_view_;
  VectorViewType mu_view_;
  ViewVectorDiscreteFunctionType u_;
  ViewDiscreteFunctionType p_;
  std::vector<ViewVectorDiscreteFunctionType> P_;
  std::vector<ViewVectorDiscreteFunctionType> Pnat_;
  std::vector<ViewDiscreteFunctionType> phi_;
  std::vector<ViewDiscreteFunctionType> phinat_;
  std::vector<ViewDiscreteFunctionType> mu_;
  const double Re_;
  const double Fa_inv_;
  const double xi_;
  const double kappa_;
  const double c_1_;
  const double Pa_;
  const double beta_;
  const double gamma_;
  const double Be_;
  const double Ca_;
  const double epsilon_;
  const double In_;
  const double vol_domain_;
  MatrixType S_stokes_;
  MatrixViewType A_stokes_;
  MatrixViewType B_stokes_;
  MatrixViewType BT_stokes_;
  MatrixType M_p_stokes_;
  EigenVectorType stokes_rhs_vector_;
  EigenVectorViewType stokes_f_vector_;
  EigenVectorViewType stokes_g_vector_;
  VectorType p_basis_integrated_vector_;
  const XT::Grid::AllDirichletBoundaryInfo<PI> boundary_info_;
  DirichletConstraints<PI, SpaceInterface<PGV, d, 1, R>> u_dirichlet_constraints_;
  DirichletConstraints<PI, SpaceInterface<PGV, 1, 1, R>> phi_dirichlet_constraints_;
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, d>> A_stokes_operator_;
  XT::LA::SparsityPatternDefault ofield_submatrix_pattern_;
  MatrixType S_ofield_;
  MatrixType M_ofield_;
  MatrixType C_ofield_elliptic_part_;
  MatrixType C_ofield_linear_part_;
  mutable MatrixViewType S_ofield_00_;
  mutable MatrixViewType S_ofield_01_;
  mutable MatrixViewType S_ofield_10_;
  mutable MatrixViewType S_ofield_11_;
  mutable std::shared_ptr<MatrixOperator<MatrixViewType, PGV, d>> S_ofield_00_operator_;
  mutable std::shared_ptr<MatrixOperator<MatrixViewType, PGV, d>> S_ofield_10_operator_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> M_ofield_operator_;
  VectorType ofield_rhs_vector_;
  mutable std::vector<EigenVectorType> ofield_old_result_;
  VectorViewType ofield_f_vector_;
  VectorViewType ofield_g_vector_;
  VectorType ofield_residual_;
  VectorType ofield_x_n_;
  VectorType ofield_update_;
  VectorType ofield_candidate_;
  mutable ColMajorBackendType S_colmajor_;
  std::shared_ptr<StokesSolverType> stokes_solver_;
  mutable std::shared_ptr<SolverType> ofield_solver_;
  mutable std::shared_ptr<SolverType> pfield_solver_;
  const bool linearize_;
  const XT::LA::SparsityPatternDefault pfield_submatrix_pattern_;
  MatrixType S_pfield_;
  MatrixType M_pfield_;
  MatrixType E_pfield_;
  MatrixType A_pfield_linear_part_;
  MatrixType A_pfield_nonlinear_part_;
  MatrixType J_pfield_linear_part_;
  MatrixViewType S_pfield_00_;
  MatrixViewType S_pfield_01_;
  MatrixViewType S_pfield_10_;
  MatrixViewType S_pfield_11_;
  MatrixViewType S_pfield_12_;
  MatrixViewType S_pfield_20_;
  MatrixViewType S_pfield_22_;
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, 1>> S_pfield_00_operator_;
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, 1>> S_pfield_10_operator_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> A_pfield_nonlinear_part_operator_;
  VectorType pfield_rhs_vector_;
  std::vector<EigenVectorType> pfield_old_result_;
  XT::LA::VectorView<VectorType> pfield_f_vector_;
  XT::LA::VectorView<VectorType> pfield_g_vector_;
  XT::LA::VectorView<VectorType> pfield_h_vector_;
  mutable VectorDiscreteFunctionType u_discr_func_;
  mutable std::vector<VectorDiscreteFunctionType> P_discr_func_;
  mutable std::vector<VectorDiscreteFunctionType> Pnat_discr_func_;
  mutable std::vector<DiscreteFunctionType> phi_discr_func_;
  mutable std::vector<DiscreteFunctionType> phinat_discr_func_;
  mutable std::vector<DiscreteFunctionType> mu_discr_func_;
  mutable std::shared_ptr<StokesPerThread> u_discr_func_local_;
  mutable std::shared_ptr<OfieldPerThread> P_discr_func_local_;
  mutable std::shared_ptr<OfieldPerThread> Pnat_discr_func_local_;
  mutable std::shared_ptr<PfieldPerThread> phi_discr_func_local_;
  mutable std::shared_ptr<PfieldPerThread> phinat_discr_func_local_;
  mutable std::shared_ptr<PfieldPerThread> mu_discr_func_local_;
};

#endif // DUNE_GDT_EXAMPLES_CELLMODEL_HH