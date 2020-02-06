// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_CELLMODEL_SOLVER_HH
#define DUNE_GDT_TEST_CELLMODEL_SOLVER_HH

#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/vector-view.hh>
#include <dune/xt/la/container/matrix-view.hh>

#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>

#include "linearoperators.hh"
#include "linearsolvers.hh"
#include "preconditioners.hh"
#include "scalarproducts.hh"

using namespace Dune::GDT;

namespace Dune {


// forward declarations
template <class X, class Y, class F>
class FGMResSolver;


struct CellModelSolver
{
  using G = ALU_2D_SIMPLEX_CONFORMING;
  // using G = YASP_2D_EQUIDISTANT_OFFSET;
  static const constexpr size_t d = G::dimension;
  using GV = typename G::LeafGridView;
  using PGV = XT::Grid::PeriodicGridView<GV>;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using PI = XT::Grid::extract_intersection_t<PGV>;
  using MatrixType = XT::LA::EigenRowMajorSparseMatrix<double>;
  using VectorType = XT::LA::CommonDenseVector<double>;
  using EigenVectorType = XT::LA::EigenDenseVector<double>;
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
  using LUSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using OfieldDirectSolverType = ::Eigen::SparseLU<ColMajorBackendType>;
  using PhiDirichletConstraintsType = DirichletConstraints<PI, SpaceInterface<PGV, 1, 1, R>>;
  using OfieldSchurSolverType = Dune::RestartedGMResSolver<EigenVectorType>;
  using PerThreadVectorLocalFunc = XT::Common::PerThreadValue<std::unique_ptr<VectorLocalDiscreteFunctionType>>;
  using PerThreadScalarLocalFuncs = XT::Common::PerThreadValue<std::vector<std::unique_ptr<LocalDiscreteFunctionType>>>;
  using PerThreadVectorLocalFuncs =
      XT::Common::PerThreadValue<std::vector<std::unique_ptr<VectorLocalDiscreteFunctionType>>>;

  CellModelSolver(
      const std::string testcase = "single_cell",
      const double t_end = 1.,
      const unsigned int num_elements_x = 50,
      const unsigned int num_elements_y = 50,
      const bool use_tbb = true,
      const double Be = 0.3, // bending capillary number, ratio of viscous forces to bending forces
      const double Ca = 0.1, // capillary number, ratio of viscous forces to surface tension forces
      const double Pa = 1, // polarization elasticity number
      const double Re = 5e-13, // Reynolds number
      const double Fa = 1., // active force number
      const double xi = 1.1, // alignment of P with the flow, > 0 for rod-like cells and < 0 for oblate ones
      const double kappa = 1.65, // eta_rot/eta, scaling factor between rotational and dynamic viscosity
      const double c_1 = 5., // double well shape parameter
      const double beta = 0., // alignment of P with the boundary of cell
      const double gamma = 0.025, // phase field mobility coefficient
      const double epsilon = 0.21, // phase field parameter
      const double In = 1., // interaction parameter
      const CellModelLinearSolverType pfield_solver_type = CellModelLinearSolverType::schur_fgmres_gmres,
      const CellModelMassMatrixSolverType pfield_mass_matrix_solver_type = CellModelMassMatrixSolverType::sparse_lu,
      const CellModelLinearSolverType ofield_solver_type = CellModelLinearSolverType::schur_fgmres_gmres,
      const CellModelMassMatrixSolverType ofield_mass_matrix_solver_type = CellModelMassMatrixSolverType::sparse_lu,
      const double outer_reduction = 1e-10,
      const int outer_restart = 100,
      const int outer_verbose = 0,
      const double inner_reduction = 1e-3,
      const int inner_maxit = 10,
      const int inner_verbose = 0,
      const bool linearize = false,
      const double pol_order = 1);

  size_t num_cells() const;

  bool linear() const;

  bool finished() const;

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
  std::vector<std::vector<VectorType>> solve(const double dt,
                                             const bool write,
                                             const double write_step,
                                             const std::string filename = "cellmodel",
                                             const bool subsampling = true);

  // Like solve, but only computes and returns the next n timesteps
  std::vector<std::vector<VectorType>> next_n_timesteps(const size_t n, const double dt);

  //******************************************************************************************************************
  //********************************* Product operators (mass matrix application) ************************************
  //******************************************************************************************************************

  // applies the pfield mass matrix to phi, phinat, mu
  // To calculate the sum of the squared L2 products of phi, phinat and mu, calculate the inner product of the result
  // with vec.
  VectorType apply_pfield_product_operator(const VectorType& vec) const;

  // applies the ofield mass matrix to P, Pnat
  // To calculate the sum of the squared L2 products of P and Pnat, calculate the inner product of the result with vec.
  VectorType apply_ofield_product_operator(const VectorType& vec) const;

  VectorType apply_stokes_product_operator(const VectorType& vec) const;

  //******************************************************************************************************************
  //*****************************************  Visualization   *******************************************************
  //******************************************************************************************************************

  // Visualizes given vector as phasefield finite element vector
  void visualize_pfield(const std::string& filename, const VectorType& vec, const bool subsampling = true) const;

  // Visualizes given vector as orientation field finite element vector
  void visualize_ofield(const std::string& filename, const VectorType& vec, const bool subsampling = true) const;

  // Visualizes given vector as stokes finite element vector
  void visualize_stokes(const std::string& filename, const VectorType& vec, const bool subsampling = true) const;

  // Visualizes variables currently stored in this class.
  // If txt = true, also writes textfiles containing the values.
  void visualize(const std::string& prefix,
                 const size_t step,
                 const double t,
                 const bool subsampling = true,
                 const bool vtu = true,
                 const bool txt = false) const;

  //******************************************************************************************************************
  //*******************************  Methods to get and set variable values   ****************************************
  //******************************************************************************************************************

  // Sets stokes vector to stokes_vec
  void set_stokes_vec(const VectorType& stokes_vec);

  // Sets orientation field vector belonging to cell to pfield_vec
  void set_ofield_vec(const size_t cell, const VectorType& ofield_vec);

  // Sets phasefield vector belonging to cell to pfield_vec
  void set_pfield_vec(const size_t cell, const VectorType& pfield_vec);

  // Get stokes finite element vector
  const VectorType& stokes_vec();

  // Get orientation field finite element vector belonging to cell
  const VectorType& ofield_vec(const size_t cell);

  // Get phase field finite element vector belonging to cell
  const VectorType& pfield_vec(const size_t cell);

  //******************************************************************************************************************
  //****** Prepare methods (calculate everything that is linear for the respective operator, but depends on **********
  //****** the values of other variables, so cannot be computed once and for all in the constructor )       **********
  //******************************************************************************************************************

  void prepare_stokes_operator();

  void prepare_ofield_operator(const double dt, const size_t cell, const bool restricted = false);

  void prepare_pfield_operator(const double dt, const size_t cell, const bool restricted = false);

  void compute_restricted_ofield_dofs(const std::vector<size_t>& output_dofs, const size_t cell);

  void compute_restricted_pfield_dofs(const std::vector<size_t>& output_dofs, const size_t cell);

  //******************************************************************************************************************
  //*********************************************** Apply operators **************************************************
  //******************************************************************************************************************

  // Applies stokes operator (applies the F if Stokes equation is F(y) = 0)
  VectorType apply_stokes_operator(VectorType y, const bool /*restricted*/ = false) const;

  void assemble_nonlinear_part_of_ofield_residual(VectorType& residual, const size_t cell, const bool restricted);

  // Applies cell-th orientation field operator (applies F if the orientation field equation is F(y) = 0)
  VectorType apply_ofield_operator(const VectorType& y, const size_t cell, const bool restricted = false);

  void update_ofield_parameters(const double Pa);

  // Applies cell-th phase field operator (applies F if phase field equation is F(y) = 0)
  VectorType apply_pfield_operator(const VectorType& y, const size_t cell, const bool restricted = false);

  void update_pfield_parameters(const double Be, const double Ca, const double Pa);

  //******************************************************************************************************************
  //******************************************* Apply inverse operators **********************************************
  //******************************************************************************************************************

  // Applies inverse stokes operator (solves F(y) = 0)
  VectorType apply_inverse_stokes_operator() const;

  // Applies inverse orientation field operator (solves F(y) = 0)
  // y_guess is the initial guess for the Newton iteration
  VectorType apply_inverse_ofield_operator(const VectorType& y_guess, const size_t cell);

  // Applies inverse phase field operator (solves F(y) = 0)
  // y_guess is the initial guess for the Newton iteration
  VectorType apply_inverse_pfield_operator(const VectorType& y_guess, const size_t cell);

  //******************************************************************************************************************
  //********************************************** Apply jacobians ***************************************************
  //******************************************************************************************************************

  void set_pfield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted = false);

  // Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
  // As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
  VectorType apply_pfield_jacobian(const VectorType& source, const size_t cell, const bool restricted = false);

  VectorType apply_inverse_pfield_jacobian(const VectorType& rhs, const size_t cell);

  void set_ofield_jacobian_state(const VectorType& source, const size_t cell, const bool restricted = false);

  // Currently takes a full-dimensional vector, but only applies the rows that are in pfield_output_dofs
  // As the rows are sparse, there shouldn't be too much performance impact of applying to the whole vector
  VectorType apply_ofield_jacobian(const VectorType& source, const size_t cell, const bool restricted = false);

  VectorType apply_inverse_ofield_jacobian(const VectorType& rhs, const size_t cell);

  VectorType apply_stokes_jacobian(const VectorType& source, const bool /*restricted*/ = false);

  VectorType apply_inverse_stokes_jacobian(const VectorType& rhs);

  //******************************************************************************************************************
  //**************************** Methods to assemble rhs, residuals and jacobians ************************************
  //******************************************************************************************************************

  // Computes stokes rhs using currently stored values of variables and stores in stokes_rhs_vector_
  void assemble_stokes_rhs();

  // Computes orientation field rhs using currently stored values of variables and stores in ofield_rhs_vector_
  void assemble_ofield_rhs(const double /*dt*/, const size_t cell);

  // Computes phase field rhs using currently stored values of variables and stores in pfield_rhs_vector_
  void assemble_pfield_rhs(const double /*dt*/, const size_t cell, const bool restricted);

  // assembles linear part of orientation field jacobian and stores in S_ofield_
  void assemble_ofield_linear_jacobian(const double dt, const size_t cell);

  // assembles nonlinear part of orientation field jacobian and adds to S_ofield_
  // if assemble_ofield_linear_jacobian has been called first, S_ofield now contains the whole orientation field
  // jacobian
  void assemble_ofield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted = false);

  void assemble_C_ofield_nonlinear_part(const size_t cell, const bool restricted = false);

  // assembles linear part of phase field jacobian
  void assemble_pfield_linear_jacobian(const double /*dt*/, const size_t cell, const bool restricted);

  void set_mat_to_zero(MatrixType& mat,
                       const bool restricted,
                       const XT::LA::SparsityPatternDefault& pattern,
                       const std::vector<size_t>& rows);

  void assemble_D_pfield(const size_t cell, const bool restricted);

  // assembles nonlinear part of phase field jacobian
  void assemble_pfield_nonlinear_jacobian(const VectorType& y, const size_t cell, const bool restricted = false);

  // stores matrix with entries \int (3 phi^2 - 1) varphi_i varphi_j in M_nonlin_pfield_
  void assemble_M_nonlin_pfield(const size_t cell, const bool restricted);

  // stores nonlinear part of block G of the phase field jacobian matrix in G_pfield_nonlinear_part_
  void assemble_G_pfield(const size_t cell, const bool restricted);

  // assembles nonlinear part of phasefield residual and adds to residual
  void assemble_nonlinear_part_of_pfield_residual(VectorType& residual, const size_t cell, const bool restricted);

  //******************************************************************************************************************
  //******************************************* DEIM related methods *************************************************
  //******************************************************************************************************************

  // Dofs needed for evaluation of output_dofs provided in
  std::vector<size_t> pfield_deim_input_dofs(const size_t cell) const;

  size_t pfield_deim_input_dofs_size(const size_t cell) const;

  // Dofs needed for evaluation of output_dofs provided in
  std::vector<size_t> ofield_deim_input_dofs(const size_t cell) const;

  // private:
  //******************************************************************************************************************
  //************ The following methods all bind or evaluate the respective temporary discrete function ***************
  //******************************************************************************************************************

  void bind_u(const E& element) const;

  DomainRetType eval_u(const DomainType& x_local, const XT::Common::Parameter& param) const;

  JacobianRetType grad_u(const DomainType& x_local, const XT::Common::Parameter& param) const;

  void bind_P(const size_t cell, const E& element) const;

  DomainRetType eval_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const;

  JacobianRetType grad_P(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const;

  void bind_Pnat(const size_t cell, const E& element) const;

  DomainRetType eval_Pnat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param) const;

  void bind_phi(const size_t cell, const E& element) const;

  R eval_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param);

  DomainRetType grad_phi(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param);

  void bind_phinat(const size_t cell, const E& element) const;

  R eval_phinat(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param);

  void bind_mu(const size_t cell, const E& element) const;

  R eval_mu(const size_t cell, const DomainType& x_local, const XT::Common::Parameter& param);

  //******************************************************************************************************************
  //****************  Linear algebra operations acting only on parts of the given matrices and vectors ***************
  //******************************************************************************************************************

  // Matrix vector multiplication.
  // In the restricted case, only uses the matrix rows provided and does not touch the other entries of range.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(const MatrixType&, const VecType1&, VecType2&, const std::vector<size_t>&)>
  mv_func(const bool restricted) const
  {
    if (!restricted) {
      return [](const MatrixType& mat, const VecType1& source, VecType2& range, const std::vector<size_t>&) {
        mat.mv(source, range);
      };
    } else {
      return [](const MatrixType& mat, const VecType1& source, VecType2& range, const std::vector<size_t>& rows) {
        for (const auto& row : rows) {
          range[row] = 0.;
          for (typename MatrixType::BackendType::InnerIterator it(mat.backend(), row); it; ++it)
            range[row] += it.value() * source[it.col()];
        }
      };
    } // if (restricted)
  }

  // Multiplies given entries of vec by alpha.
  template <class VecType>
  std::function<void(VecType&, const R, const std::vector<size_t>&)> scal_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType& vec, const R alpha, const std::vector<size_t>&) { vec *= alpha; };
    } else {
      return [](VecType& vec, const R alpha, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          vec[dof] *= alpha;
      };
    } // if (restricted)
  }

  // Adds given entries of rhs to respective entries of lhs.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const VecType2&, const std::vector<size_t>&)> copy_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>&) { lhs = rhs; };
    } else {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] = rhs[dof];
      };
    } // if (restricted)
  }

  // Adds given entries of rhs to respective entries of lhs.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const VecType2&, const std::vector<size_t>&)> add_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>&) { lhs += rhs; };
    } else {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] += rhs[dof];
      };
    } // if (restricted)
  }

  // Subtracts given entries of rhs from respective entries of lhs.
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const VecType2&, const std::vector<size_t>&)> sub_func(const bool restricted) const
  {
    if (!restricted) {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>&) { lhs -= rhs; };
    } else {
      return [](VecType1& lhs, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] -= rhs[dof];
      };
    } // if (restricted)
  }

  // Computes lhs += alpha * rhs;
  template <class VecType1, class VecType2 = VecType1>
  std::function<void(VecType1&, const R, const VecType2&, const std::vector<size_t>&)>
  axpy_func(const bool restricted) const
  {
    if (!restricted) {
      return
          [](VecType1& lhs, const R alpha, const VecType2& rhs, const std::vector<size_t>&) { lhs.axpy(alpha, rhs); };
    } else {
      return [](VecType1& lhs, const R alpha, const VecType2& rhs, const std::vector<size_t>& dofs) {
        for (const auto& dof : dofs)
          lhs[dof] += rhs[dof] * alpha;
      };
    } // if (restricted)
  }

  // Copies low-dimensional vec to given entries of high-dimensional vec.
  void copy_ld_to_hd_vec(const std::vector<size_t> dofs, const VectorType& ld_vec, VectorType& hd_vec);

  //***************************************************************************************************************
  //*********************************************  Helper methods  ************************************************
  //***************************************************************************************************************

  // get lower left of computational domain from testcase name
  static XT::Common::FieldVector<R, d> get_lower_left(const std::string& testcase);

  // get upper right of computational domain from testcase name
  static XT::Common::FieldVector<R, d> get_upper_right(const std::string& testcase);

  // get directions in which domain is periodic from testcase name
  static std::string get_periodic_directions(const std::string& testcase);

  // get number of cells from testcase name
  static size_t get_num_cells(const std::string& testcase);

  // creates sparsity pattern of stokes system matrix
  static XT::LA::SparsityPatternDefault create_stokes_pattern(const SpaceInterface<PGV, d, 1, R>& u_space,
                                                              const SpaceInterface<PGV, 1, 1, R>& p_space);

  // appends entity to input_entities if one of its global_indices is in output_dofs
  void maybe_add_entity(const E& entity,
                        const DynamicVector<size_t>& global_indices,
                        const std::vector<size_t>& output_dofs,
                        std::vector<E>& input_entities,
                        const size_t subvector_size) const;

  // sets temporary orientation field discrete functions to source values
  void fill_tmp_ofield(const size_t cell, const VectorType& source, const bool restricted = false) const;

  // sets temporary phase field discrete functions to source values
  void fill_tmp_pfield(const size_t cell, const VectorType& source, const bool restricted) const;

  // error norm used in orientation field Newton iteration
  // TODO: use appropriate norm
  double ofield_residual_norm(const VectorType& residual, double l2_ref_P, double l2_ref_Pnat) const;

  // error norm used in phase field Newton iteration
  // TODO: use appropriate norm
  double
  pfield_residual_norm(const VectorType& residual, double l2_ref_phi, double l2_ref_phinat, double l2_ref_mu) const;

  R B_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param);

  R w_func(const size_t kk, const DomainType& x_local, const XT::Common::Parameter& param);

  //******************************************************************************************************************
  //*******************************************  Member variables ****************************************************
  //******************************************************************************************************************

  // Model parameters
  XT::Common::FieldVector<R, d> lower_left_;
  XT::Common::FieldVector<R, d> upper_right_;
  const double t_end_;
  double t_;
  const bool use_tbb_;
  const double Re_;
  double Fa_inv_;
  double xi_;
  double kappa_;
  double c_1_;
  double Pa_;
  double last_pfield_Pa_;
  double last_ofield_Pa_;
  double beta_;
  double gamma_;
  double Be_;
  double Ca_;
  double epsilon_;
  double In_;
  const double vol_domain_;
  const size_t num_cells_;
  const bool linearize_;
  // Grid and grid views
  XT::Grid::GridProvider<G> grid_;
  const GV nonperiodic_grid_view_;
  const PGV grid_view_;
  // Finite element function spaces
  const ContinuousLagrangeSpace<PGV, d, R> u_space_;
  const ContinuousLagrangeSpace<PGV, d, R> P_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> p_space_;
  const ContinuousLagrangeSpace<PGV, 1, R> phi_space_;
  // Size of finite element vectors
  const size_t size_u_;
  const size_t size_P_;
  const size_t size_p_;
  const size_t size_phi_;
  const size_t num_mutexes_u_;
  const size_t num_mutexes_ofield_;
  const size_t num_mutexes_pfield_;
  // Finite element vectors for phase field, orientation field and stokes system
  // There is one phase field and orientation field per cell
  VectorType stokes_vector_;
  std::vector<VectorType> ofield_vectors_;
  std::vector<VectorType> pfield_vectors_;
  // Views on parts of the system vectors corresponding to the respective variables
  VectorViewType u_view_;
  VectorViewType p_view_;
  std::vector<VectorViewType> P_view_;
  std::vector<VectorViewType> Pnat_view_;
  std::vector<VectorViewType> phi_view_;
  std::vector<VectorViewType> phinat_view_;
  std::vector<VectorViewType> mu_view_;
  // DiscreteFunctions corresponding to view vectors
  ViewVectorDiscreteFunctionType u_;
  ViewDiscreteFunctionType p_;
  std::vector<ViewVectorDiscreteFunctionType> P_;
  std::vector<ViewVectorDiscreteFunctionType> Pnat_;
  std::vector<ViewDiscreteFunctionType> phi_;
  std::vector<ViewDiscreteFunctionType> phinat_;
  std::vector<ViewDiscreteFunctionType> mu_;
  // Stokes system matrix S = (A B; B^T 0) and views on matrix blocks
  MatrixType S_stokes_;
  MatrixViewType A_stokes_;
  MatrixViewType B_stokes_;
  MatrixViewType BT_stokes_;
  // pressure mass matrix
  MatrixType M_p_stokes_;
  // Matrix operator for A_stokes_
  std::shared_ptr<MatrixOperator<MatrixViewType, PGV, d>> A_stokes_op_;
  // finite element vector rhs = (f; g) for stokes system and views on velocity and pressure parts f and g
  EigenVectorType stokes_rhs_vector_;
  EigenVectorViewType stokes_f_vector_;
  EigenVectorViewType stokes_g_vector_;
  // vector containing integrals of pressure basis functions (for normalizing such that \int p = 0)
  VectorType p_basis_integrated_vector_;
  // Dirichlet constraints
  const XT::Grid::AllDirichletBoundaryInfo<PI> boundary_info_;
  DirichletConstraints<PI, SpaceInterface<PGV, d, 1, R>> u_dirichlet_constraints_;
  PhiDirichletConstraintsType phi_dirichlet_constraints_;
  // phi is shifted by this value, i.e., instead of solving for phi, we are solving for phi + phi_shift.
  // Set to 1 by default to get 0 on the boundary instead of -1.
  R phi_shift_;
  R phinat_scale_factor_;
  // Sparsity pattern of one block of orientation field system matrix
  XT::LA::SparsityPatternDefault ofield_submatrix_pattern_;
  // Orientation field mass matrix
  MatrixType M_ofield_;
  MatrixType A_ofield_;
  // Part of C that is independent of phi and dt
  MatrixType C_ofield_elliptic_part_;
  // Whole linear part of C (elliptic part + phi-dependent part)
  MatrixType C_ofield_linear_part_;
  // Whole linear part of C (elliptic part + phi-dependent part)
  mutable MatrixType C_ofield_nonlinear_part_;
  // Linear part of ofield schur matrix M/dt + A - 1/kappa C
  MatrixType S_schur_ofield_linear_part_;
  // Matrix operators for orientation field matrices
  std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> M_ofield_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> A_ofield_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> C_ofield_linear_part_op_;
  mutable std::shared_ptr<MatrixOperator<MatrixType, PGV, d>> C_ofield_nonlinear_part_op_;
  OfieldMatrixLinearPartOperator<VectorType, MatrixType, CellModelSolver> ofield_jac_linear_op_;
  OfieldLinearSolver ofield_solver_;
  // finite element vector rhs = (f; g) for ofield system and views on P and Pnat parts f and g
  VectorType ofield_rhs_vector_;
  VectorViewType ofield_f_vector_;
  VectorViewType ofield_g_vector_;
  // Linear solvers and linear operators needed for solvers
  mutable ColMajorBackendType S_colmajor_;
  std::shared_ptr<LUSolverType> stokes_solver_;
  VectorType ofield_tmp_vec_;
  VectorType ofield_tmp_vec2_;
  // Indices for restricted operator in DEIM context
  std::vector<std::vector<size_t>> ofield_deim_input_dofs_;
  std::vector<size_t> Pnat_deim_input_dofs_begin_;
  // output dofs that were computed by the DEIM algorithm
  std::vector<std::shared_ptr<std::vector<size_t>>> ofield_deim_output_dofs_;
  std::vector<std::vector<size_t>> ofield_deim_unique_output_dofs_;
  // output dofs for the respective variables, shifted to range [0, size_phi)
  std::vector<std::vector<size_t>> P_deim_output_dofs_;
  std::vector<std::vector<size_t>> Pnat_deim_output_dofs_;
  // Entities that we have to walk over to calculate values at output dofs
  std::vector<std::vector<E>> ofield_deim_entities_;
  // DiscreteFunctions and vectors to be used in (nonlinear) calculations where a source vector is provided
  // Sparsity pattern of one block of phase field system matrix
  const XT::LA::SparsityPatternDefault pfield_submatrix_pattern_;
  // Phase field system matrix S = (M/dt+D E 0; G H J; A 0 C)
  MatrixType M_pfield_;
  MatrixType D_pfield_;
  MatrixType M_ell_pfield_;
  MatrixType M_nonlin_pfield_;
  MatrixType G_pfield_;
  MatrixType A_boundary_pfield_;
  mutable ColMajorBackendType M_pfield_colmajor_;
  mutable ColMajorBackendType S_pfield_colmajor_;
  // Pfield solvers and linear operators
  // Matrix operators for phasefield matrices
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> D_pfield_op_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> G_pfield_op_;
  std::shared_ptr<MatrixOperator<MatrixType, PGV, 1>> M_nonlin_pfield_op_;
  PfieldMatrixLinearPartOperator<VectorType, MatrixType, PhiDirichletConstraintsType, CellModelSolver>
      pfield_jac_linear_op_;
  PfieldLinearSolver pfield_solver_;
  // Phase field rhs vector (r0 r1 r2)
  VectorType pfield_rhs_vector_;
  VectorViewType pfield_r0_vector_;
  VectorViewType pfield_r1_vector_;
  VectorViewType pfield_r2_vector_;
  mutable EigenVectorType phi_tmp_eigen_;
  mutable EigenVectorType phi_tmp_eigen2_;
  mutable EigenVectorType phi_tmp_eigen3_;
  mutable EigenVectorType phi_tmp_eigen4_;
  // Indices for restricted operator in DEIM context
  std::vector<std::vector<size_t>> pfield_deim_input_dofs_;
  // phinat_deim_input_dofs_begin_[cell] contains index of first phinat input dof in pfield_deim_input_dofs_[cell]
  // vector
  std::vector<size_t> phinat_deim_input_dofs_begin_;
  std::vector<size_t> mu_deim_input_dofs_begin_;
  // output dofs that were computed by the DEIM algorithm
  std::vector<std::shared_ptr<std::vector<size_t>>> pfield_deim_output_dofs_;
  std::vector<std::vector<size_t>> pfield_deim_unique_output_dofs_;
  // output dofs for the respective variables, shifted to range [0, size_phi)
  std::vector<std::vector<size_t>> phi_deim_output_dofs_;
  std::vector<std::vector<size_t>> phinat_deim_output_dofs_;
  std::vector<std::vector<size_t>> mu_deim_output_dofs_;
  // Indices that are both in phinat and mu output dofs
  std::vector<std::vector<size_t>> both_mu_and_phi_deim_output_dofs_;
  // Entities that we have to walk over to calculate values at output dofs
  std::vector<std::vector<E>> pfield_deim_entities_;
  // DiscreteFunctions and vectors to be used in (nonlinear) calculations where a source vector is provided
  VectorType pfield_tmp_vec_;
  mutable EigenVectorType pfield_tmp_eigen_;
  VectorType pfield_tmp_vec2_;
  VectorType phi_tmp_vec_;
  VectorType phi_tmp_vec2_;
  VectorType u_tmp_vec_;
  VectorType P_tmp_vec_;
  mutable VectorDiscreteFunctionType u_tmp_;
  mutable std::vector<VectorDiscreteFunctionType> P_tmp_;
  mutable std::vector<VectorDiscreteFunctionType> Pnat_tmp_;
  mutable std::vector<DiscreteFunctionType> phi_tmp_;
  mutable std::vector<DiscreteFunctionType> phinat_tmp_;
  mutable std::vector<DiscreteFunctionType> mu_tmp_;
  mutable std::shared_ptr<PerThreadVectorLocalFunc> u_tmp_local_;
  mutable std::shared_ptr<PerThreadVectorLocalFuncs> P_tmp_local_;
  mutable std::shared_ptr<PerThreadVectorLocalFuncs> Pnat_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> phi_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> phinat_tmp_local_;
  mutable std::shared_ptr<PerThreadScalarLocalFuncs> mu_tmp_local_;
  double dt_;
};


} // namespace Dune

#endif // DUNE_GDT_TEST_CELLMODEL_SOLVER_HH
