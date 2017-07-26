// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_IMPLICIT_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_IMPLICIT_RUNGEKUTTA_HH

#include <utility>


#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/pattern.hh>

#include <dune/gdt/assembler/codim0-matrix-datahandle.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "interface.hh"
#include "explicit-rungekutta.hh"


namespace Dune {
namespace GDT {


namespace internal {


// backward Euler
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::implicit_euler>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }
};

// implicit midpoint
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::implicit_midpoint>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0.5]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5]");
  }
};

// Trapezoidal rule
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::trapezoidal_rule>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0; 0.5 0.5]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5 0.5]");
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0 1]");
  }
};


} // namespace internal


/** \brief Time stepper using diagonally implicit Runge Kutta methods
 *
 * Timestepper using diagonally implicit Runge Kutta methods to solve equations of the form u_t = r * L(u, t) where u is
 * a
 * discrete function, L an operator acting on u and r a scalar factor (e.g. -1).
 * The specific Runge Kutta method can be chosen as the third template argument. If your desired Runge Kutta method is
 * not contained in DiagonallyImplicitRungeKuttaMethods, choose TimeStepperMethods::diagonally_implicit_other and supply
 * a
 * DynamicMatrix< RangeFieldType > A and vectors (DynamicVector< RangeFieldType >) b and c in the constructor. Here, A,
 * b and c form the butcher tableau (see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods, A is
 * composed of the coefficients a_{ij}, b of b_j and c of c_j). The default is a backward euler method.
 *
 * Notation (the factor r is omitted for notational simplicity): For an s-stage method,
 * \mathbf{u}^{n+1} = \mathbf{u}^n + dt \sum_{i=0}^{s-1} b_i \mathbf{k}_i
 * \mathbf{k}_i = L(\mathbf{u}_i, t_i)
 * t_i = t^n + dt c_i
 * \mathbf{u}_i = \mathbf{u}^n + dt \sum_{j=0}^{s-1} a_{ij} \mathbf{k}_j,
 * To solve for the \mathbf{u}_i, a Newton method with Armijo backtracking line search is used at each stage (if
 * a_{ii} \neq 0). The Newton update for each stage takes the form
 * N(\mathbf{u}_l) \mathbf{d}_l := (I - dt a_{ii} J_L(\mathbf{u}_l, t_i)) \mathbf{d}_l = -\mathbf{u}_l
 * + \mathbf{u}_{expl, i} + dt a_{ii} L(\mathbf{u}_l, t_i) =: -\mathbf{res}(\mathbf{u}_l),
 * where \mathbf{u}_l is the current iterate,
 * \mathbf{u}_{expl, i} = \mathbf{u}^n + dt \sum_{j=0}^{j-1} b_j \mathbf{k}_j
 * the explicit part of \mathbf{u}_i, J_L the jacobian of the operator L with respect to \mathbf{u} and
 * \mathbf{u}_{l+1} = \mathbf{u}_l + \alpha \mathbf{d}_l. The factor \alpha is found by an Armijo backtracking
 * line search such that
 * \mathbf{res}(\mathbf{u}_l + \alpha * \mathbf{d}_l) < \mathbf{res}(\mathbf{u}^n) + \beta D\mathbf{res} \alpha
 * \mathbf{d}_l,
 * where D\mathbf{res} = -N is the jacobian of the residual.
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values
 */
template <class OperatorImp,
          class DiscreteFunctionImp,
          TimeStepperMethods method = TimeStepperMethods::implicit_euler,
          XT::LA::Backends container_backend = XT::LA::default_sparse_backend>
class DiagonallyImplicitRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp> BaseType;
  typedef typename internal::ButcherArrayProvider<typename BaseType::RangeFieldType, method> ButcherArrayProviderType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;
  using typename BaseType::DataHandleType;

  typedef OperatorImp OperatorType;
  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::DynamicVector<RangeFieldType> VectorType;
  typedef typename XT::LA::Container<RangeFieldType, container_backend>::MatrixType SolverMatrixType;
  typedef typename XT::LA::Container<RangeFieldType, container_backend>::VectorType SolverVectorType;
  typedef typename DiscreteFunctionType::SpaceType::CommunicatorType CommunicatorType;
  typedef typename XT::LA::Solver<SolverMatrixType, CommunicatorType> SolverType;
  typedef typename Dune::GDT::MatrixDataHandle<SolverMatrixType, typename DiscreteFunctionType::SpaceType>
      SolverMatrixDataHandleType;

  using BaseType::current_solution;
  using BaseType::current_time;

  /**
   * \brief Constructor for RungeKutta time stepper
   * \param op Operator L
   * \param initial_values Discrete function containing initial values for u at time t_0.
   * \param r Scalar factor (see above, default is 1)
   * \param t_0 Initial time (default is 0)
   * \param solver_type String to choose solver_type, for possible values see SolverType::types().
   * Default is the empty string, which chooses the default solver type.
   * \param A Coefficient matrix (only provide if you use TimeStepperMethods::diagonally_implicit_other)
   * \param b Coefficient vector (only provide if you use TimeStepperMethods::diagonally_implicit_other)
   * \param c Coefficients for time steps (only provide if you use TimeStepperMethods::diagonally_implicit_other)
   */
  DiagonallyImplicitRungeKuttaTimeStepper(const OperatorType& op,
                                          const DiscreteFunctionType& initial_values,
                                          const RangeFieldType r = 1.0,
                                          const RangeFieldType t_0 = 0.0,
                                          const RangeFieldType beta = 1e-4,
                                          const std::string solver_type = "",
                                          const MatrixType& A = ButcherArrayProviderType::A(),
                                          const VectorType& b = ButcherArrayProviderType::b(),
                                          const VectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , op_(op)
    , r_(r)
    , u_i_(BaseType::current_solution())
    , u_i_expl_(u_i_)
    , d_(u_i_)
    , N_times_d_(u_i_)
    , res_(u_i_)
    , new_res_(u_i_)
    , u_i_plus_alpha_d_(u_i_)
    , solver_type_(solver_type)
    , beta_(beta)
    , newton_matrix_(u_i_.vector().size(), u_i_.vector().size(), this->newton_matrix_pattern(u_i_.space()))
    , solver_(newton_matrix_, u_i_.space().communicator())
    , A_(A)
    , b_(b)
    , c_(c)
    , num_stages_(A_.rows())
  {
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_.size() == A_.rows());
    assert(c_.size() == A_.rows());
    bool lower_triangular = true;
    bool diagonally_implicit = true;
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        if (Dune::XT::Common::FloatCmp::ne(A_[ii][jj], 0.0))
          lower_triangular = false;
        if (jj > ii && Dune::XT::Common::FloatCmp::ne(A_[ii][jj], 0.0))
          diagonally_implicit = false;
      }
    }
    if (lower_triangular)
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
                 "Matrix A is strictly lower triangular, use ExplicitRungeKuttaTimeStepper!");
    if (!diagonally_implicit)
      DUNE_THROW(Dune::NotImplemented, "Matrix A has entries above the diagonal, this is not supported yet!");
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(u_i_);
    }
  } // constructor

  RangeFieldType calculate_l2_norm(const DiscreteFunctionType& disc_func)
  {
    RangeFieldType norm = 0;
    const auto& vector = disc_func.vector();
    const auto& mapper = disc_func.space().mapper();
    for (const auto& entity : Dune::elements(disc_func.space().grid_layer(), Dune::Partitions::interiorBorder))
      for (const auto& index : mapper.globalIndices(entity))
        norm += std::pow(vector[index], 2);
    disc_func.space().grid_layer().comm().sum(norm);
    norm = std::sqrt(norm);
    return norm;
  }

  virtual RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    // calculate stages
    SolverMatrixDataHandleType newton_matrix_handle(newton_matrix_, u_n.space());
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      // calculate u_expl
      if (ii == 0)
        u_i_expl_.vector() = u_n.vector();
      else
        u_i_expl_.vector() += stages_k_[ii - 1].vector() * (actual_dt * r_ * (A_[ii][ii - 1]));
      DataHandleType stages_k_ii_handle(stages_k_[ii]);
      // choose explicit or implicit treatment
      if (XT::Common::FloatCmp::eq(A_[ii][ii], 0.0)) {
        // stage is explicit
        stages_k_[ii].vector() *= RangeFieldType(0);
        op_.apply(
            u_i_expl_, stages_k_[ii], XT::Common::Parameter({{"t", {t + actual_dt * c_[ii]}}, {"dt", {actual_dt}}}));
        stages_k_[ii].space().grid_layer().template communicate<DataHandleType>(
            stages_k_ii_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
      } else {
        // stage is implicit, Newton method with Armijo backtracking
        // initial guess for u_i: u_i_expl_ + dt a_{ii} k_{i-1}
        u_i_.vector() =
            (ii == 0) ? u_n.vector() : u_i_expl_.vector() + stages_k_[ii - 1].vector() * (actual_dt * r_ * A_[ii][ii]);
        // calculate residual with current iterate
        stages_k_[ii].vector() *= RangeFieldType(0);
        op_.apply(u_i_, stages_k_[ii], XT::Common::Parameter({{"t", {t + actual_dt * c_[ii]}}, {"dt", {actual_dt}}}));
        stages_k_[ii].space().grid_layer().template communicate<DataHandleType>(
            stages_k_ii_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
        res_.vector() = u_i_expl_.vector() - u_i_.vector() + stages_k_[ii].vector() * (actual_dt * r_ * A_[ii][ii]);
        // calculate norm of residual
        RangeFieldType current_error = calculate_l2_norm(res_);
        RangeFieldType first_error = current_error;
        const double abs_tol = 1e-12;
        const double rel_reduction = 1e-8;
        // Newton loop
        while (current_error > abs_tol && current_error > first_error * rel_reduction) {
          // assemble Newton matrix N = I - dt * a_{ii} * J_L
          newton_matrix_ *= 0;
          op_.assemble_jacobian(newton_matrix_, u_i_, t + actual_dt * c_[ii]);
          //          u_i_.space().grid_layer().template communicate<SolverMatrixDataHandleType>(
          //              newton_matrix_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
          newton_matrix_ *= -actual_dt * r_ * A_[ii][ii];
          for (size_t kk = 0; kk < newton_matrix_.rows(); ++kk)
            newton_matrix_.add_to_entry(kk, kk, 1.);
          // solve N d = res
          d_.vector() *= 0;
          solver_.apply(res_.vector(), d_.vector(), solver_type_);

          // Armijo backtracking line search
          RangeFieldType alpha = 2;
          RangeFieldType error_threshold = abs_tol;
          while (current_error > error_threshold) {
            // decrease alpha and calculate current u_i
            alpha /= 2;
            // calculate error threshold with current alpha
            N_times_d_.vector() *= 0;
            newton_matrix_.mv(d_.vector(), N_times_d_.vector());
            N_times_d_.vector() = res_.vector() - (N_times_d_.vector() * alpha * beta_);
            error_threshold = calculate_l2_norm(N_times_d_);
            // calculate residual with current alpha
            u_i_plus_alpha_d_.vector() = u_i_.vector() + d_.vector() * alpha;
            stages_k_[ii].vector() *= RangeFieldType(0);
            op_.apply(u_i_plus_alpha_d_,
                      stages_k_[ii],
                      XT::Common::Parameter({{"t", {t + actual_dt * c_[ii]}}, {"dt", {actual_dt}}}));
            stages_k_[ii].space().grid_layer().template communicate<DataHandleType>(
                stages_k_ii_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
            new_res_.vector() = u_i_expl_.vector() - u_i_plus_alpha_d_.vector()
                                + stages_k_[ii].vector() * (actual_dt * r_ * A_[ii][ii]);
            current_error = calculate_l2_norm(new_res_);
          } // Armijo backtracking
          u_i_.vector() = u_i_plus_alpha_d_.vector();
          res_.vector() = new_res_.vector();
        } // Newton
      } // if (explicit) else (implicit)
    } // loop over stages

    // calculate value of u at next time step
    for (size_t ii = 0; ii < num_stages_; ++ii)
      u_n.vector() += stages_k_[ii].vector() * (r_ * actual_dt * b_[ii]);

    // augment time
    t += actual_dt;

    return dt;
  } // ... step(...)

private:
  static XT::LA::SparsityPatternDefault newton_matrix_pattern(const typename DiscreteFunctionType::SpaceType& space)
  {
    XT::LA::SparsityPatternDefault jacobian_pattern = space.compute_pattern();
    for (size_t ii = 0; ii < jacobian_pattern.size(); ++ii)
      jacobian_pattern.insert(ii, ii);
    return jacobian_pattern;
  }

  const OperatorType& op_;
  const RangeFieldType r_;
  DiscreteFunctionType u_i_;
  DiscreteFunctionType u_i_expl_;
  DiscreteFunctionType d_;
  DiscreteFunctionType N_times_d_;
  DiscreteFunctionType res_;
  DiscreteFunctionType new_res_;
  DiscreteFunctionType u_i_plus_alpha_d_;
  const std::string solver_type_;
  const RangeFieldType beta_;
  SolverMatrixType newton_matrix_;
  SolverType solver_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector<DiscreteFunctionType> stages_k_;
  const size_t num_stages_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_IMPLICIT_RUNGEKUTTA_HH
