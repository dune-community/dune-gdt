// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_PLAYGROUND_TIMESTEPPER_ROSENBROCK_HH
#define DUNE_GDT_PLAYGROUND_TIMESTEPPER_ROSENBROCK_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include <dune/gdt/timestepper/interface.hh>


namespace Dune {
namespace GDT {


enum class RosenbrockTimeStepperMethods
{
  GRK4A,
  GRK4T,
  other
};


namespace internal {


// unspecialized
template <class RangeFieldType, class TimeFieldType,
          RosenbrockTimeStepperMethods method = RosenbrockTimeStepperMethods::other>
struct RosenbrockButcherArrayProvider
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in RosenbrockTimeStepper's constructor for this method!");
    return Dune::DynamicMatrix<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in RosenbrockTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in RosenbrockTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in RosenbrockTimeStepper's constructor for this method!");
    return Dune::DynamicVector<TimeFieldType>();
  }

  static Dune::DynamicMatrix<RangeFieldType> Gamma()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in RosenbrockTimeStepper's constructor for this method!");
    return Dune::DynamicVector<TimeFieldType>();
  }
};

// GRK4A, see Kaps, Rentrop (1979), "Generalized Runge-Kutta methods of order four with stepsize control for stiff
// ordinary differential equations"
template <class RangeFieldType, class TimeFieldType>
struct RosenbrockButcherArrayProvider<RangeFieldType, TimeFieldType, RosenbrockTimeStepperMethods::GRK4A>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(std::string("[0 0 0 0;") + " 0.438 0 0 0;"
                                                                             + " 0.796920457938 0.0730795420615 0 0;"
                                                                             + " 0.796920457938 0.0730795420615 0 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[0.199293275701 0.482645235674 0.0680614886256 0.25]");
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[0.346325833758  0.285693175712 0.367980990530 0]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0 0.438 0.87 0.87]");
  }

  static Dune::DynamicMatrix<RangeFieldType> Gamma()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(
        std::string("[0.395  0 0 0;") + " -0.767672395484 0.395  0 0;" + " -0.851675323742  0.522967289188 0.395  0;"
        + " 0.288463109545 0.0880214273381 -0.337389840627 0.395]");
  }
};

// GRK4T
template <class RangeFieldType, class TimeFieldType>
struct RosenbrockButcherArrayProvider<RangeFieldType, TimeFieldType, RosenbrockTimeStepperMethods::GRK4T>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(std::string("[0 0 0 0;") + " 0.462 0 0 0;"
                                                                             + " -0.0815668168327 0.961775150166 0 0;"
                                                                             + " -0.0815668168327 0.961775150166 0 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<double>>(
        "[0.217487371653 0.486229037990 0 0.296283590357]");
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[-0.717088504499 1.77617912176 -0.0590906172617 0]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0 0.462 0.88020833333 0.88020833333]");
  }

  static Dune::DynamicMatrix<RangeFieldType> Gamma()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(
        std::string("[0.231 0 0 0;") + " -0.270629667752 0.231 0 0;" + " 0.311254483294 0.00852445628482 0.231 0;"
        + " 0.282816832044 -0.457959483281 -0.111208333333 0.231]");
  }
}; // GRK4T


} // namespace internal


/** \brief Time stepper using Rosenbrock-type methods
 *
 * Timestepper using Rosenbrock-type methods methods to solve equations of the form u_t = r * L(u) where u is a
 * discrete function, L an operator acting on u and r a scalar factor (e.g. -1).
 * The specific Rosenbrock-type method can be chosen as the third template argument. If your desired Rosenbrock-type
 * method is not contained in Dune::GDT::TimeStepper::RosenbrockTimeStepperMethods, choose
 * RosenbrockTimeStepperMethods::other and
 * supply matrices A, Gamma (DynamicMatrix< RangeFieldType >) and vectors b_1, b_2 (DynamicVector< RangeFieldType >) and
 * c (DynamicVector< TimeFieldType >) in the constructor. Here, A, Gamma, b_1, b_2 and c form the extended butcher
 * tableau (see https://de.wikipedia.org/wiki/Rosenbrock-Wanner-Verfahren, b_1 and b_2 are the same as for adaptive
 * Runge-Kutta schemes). The default is the GRK4T method.
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values and solution at a fixed time
 * \tparam SolverImp Type of solver used for inversion of matrix in each time step.
 * \tparam TimeFieldImp Type used for representation of time (default is double)
 * \tparam method Rosenbrock-type method that is used (default is RosenbrockTimeStepperMethods::GRK4T)
 *
 * \todo Implement concept of jacobian/time derivative of operator and finish implementation of this method.
 */
template <class OperatorImp, class DiscreteFunctionImp, class SolverImp, class TimeFieldImp = double,
          RosenbrockTimeStepperMethods method = RosenbrockTimeStepperMethods::GRK4T>
class RosenbrockTimeStepper : public TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp> BaseType;
  typedef typename internal::RosenbrockButcherArrayProvider<typename BaseType::RangeFieldType, TimeFieldImp, method>
      ButcherArrayProviderType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;

  typedef OperatorImp OperatorType;
  typedef SolverImp SolverType;
  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::DynamicVector<RangeFieldType> VectorType;
  typedef typename Dune::DynamicVector<TimeFieldType> TimeVectorType;

  using BaseType::current_solution;
  using BaseType::current_time;

  /**
   * \brief Constructor for RosenbrockTimeStepper time stepper
   * \param op Operator L
   * \param initial_values Discrete function containing initial values for u at time t_0.
   * \param r Scalar factor (see above, default is 1)
   * \param t_0 Initial time (default is 0)
   * \param tol Error tolerance for the adaptive scheme (default is 1e-4)
   * \param scale_factor_min Minimum allowed factor for time step scaling (default is 0.2).
   * \param scale_factor_max Maximum allowed factor for time step scaling (default is 5).
   * \param A Coefficient matrix (only provide if you use RosenbrockTimeStepperMethods::other)
   * \param b_1 First set of coefficients (only provide if you use RosenbrockTimeStepperMethods::other)
   * \param b_2 Second set of coefficients (only provide if you use RosenbrockTimeStepperMethods::other)
   * \param c Coefficients for time steps (only provide if you use RosenbrockTimeStepperMethods::other)
   * \param Gamma Coefficient matrix (only provide if you use RosenbrockTimeStepperMethods::other)
   */
  RosenbrockTimeStepper(const OperatorType& op, const DiscreteFunctionType& initial_values,
                        const RangeFieldType r = 1.0, const double t_0 = 0.0, const RangeFieldType tol = 1e-4,
                        const TimeFieldType scale_factor_min = 0.2, const TimeFieldType scale_factor_max = 5,
                        const MatrixType& A     = ButcherArrayProviderType::A(),
                        const VectorType& b_1   = ButcherArrayProviderType::b_1(),
                        const VectorType& b_2   = ButcherArrayProviderType::b_2(),
                        const TimeVectorType& c = ButcherArrayProviderType::c(),
                        const MatrixType& Gamma = ButcherArrayProviderType::Gamma())
    : BaseType(t_0, initial_values)
    , op_(op)
    , r_(r)
    , tol_(tol)
    , scale_factor_min_(scale_factor_min)
    , scale_factor_max_(scale_factor_max)
    , u_tmp_(BaseType::current_solution())
    , A_(A)
    , m_1_(b_1)
    , m_2_(b_2)
    , c_(c)
    , d_(Gamma.rows())
    , Gamma_(Gamma)
    , m_diff_(b_1)
    , num_stages_(A_.rows())
    , gamma_ii_equal_for_all_i_(true)
  {
    assert(Dune::XT::Common::FloatCmp::gt(tol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::le(scale_factor_min_, 1.0));
    assert(Dune::XT::Common::FloatCmp::ge(scale_factor_max_, 1.0));
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(Gamma_.rows() == Gamma_.cols() && "Gamma has to be a square matrix");
    assert(Gamma_.rows() == A_.rows() && "Sizes of A and Gamma have to match!");
    assert(b_1.size() == A_.rows());
    assert(b_2.size() == A_.rows());
    assert(c_.size() == A_.rows());
#ifndef NDEBUG
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      TimeFieldType c_calculated = 0;
      for (size_t jj = 0; jj < ii; ++jj)
        c_calculated += A_[ii][jj];
      assert(Dune::XT::Common::FloatCmp::eq(c_calculated, c_[ii]));
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        assert(Dune::XT::Common::FloatCmp::eq(A_[ii][jj], 0.0)
               && "A has to be a lower triangular matrix with 0 on the main diagonal!");
        if (jj == ii)
          assert(Dune::XT::Common::FloatCmp::ne(Gamma_[ii][jj], 0.0)
                 && "The diagonal entries of Gamma must not vanish!");
        else
          assert(Dune::XT::Common::FloatCmp::eq(Gamma_[ii][jj], 0.0) && "Gamma has to be a lower triangular matrix!");
      }
    }
#endif // NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_intermediate_stages_.emplace_back(current_solution());
    }
    // transform variables for faster calculations
    // C_ = diag(1/gamma_11, ... , 1/gamma_nn) - Gamma^(-1), A_ = A*Gamma^(-1), m = (b_1,...,b_n)*Gamma^(-1), see e.g.
    // Hairer, Wanner (1996), Solving ordinary differential equations II: Stiff and differential-algebraic problems,
    // pp 119ff
    auto Gamma_inv = Gamma_;
    Gamma_inv.invert();
    C_ = Gamma_inv;
    C_ *= -1.0;
    for (size_t ii = 0; ii < C_.rows(); ++ii)
      C_[ii][ii] += 1.0 / (Gamma_[ii][ii]);
    A_.rightmultiply(Gamma_inv);
    Gamma_inv.mtv(b_1, m_1_);
    Gamma_inv.mtv(b_2, m_2_);

    auto gamma = Gamma[0][0];
    for (size_t ii = 0; ii < Gamma.rows(); ++ii) {
      if (Dune::XT::Common::FloatCmp::ne(gamma, Gamma_[ii]))
        gamma_ii_equal_for_all_i_ = false;
      d_[ii]                      = 0.0;
      for (size_t jj = 0; jj <= ii; ++jj)
        d_[ii] += Gamma[ii][jj];
    }
  } // constructor RosenbrockTimeStepper

  TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt)
  {
    TimeFieldType actual_dt              = std::min(dt, max_dt);
    RangeFieldType mixed_error           = std::numeric_limits<RangeFieldType>::max();
    TimeFieldType time_step_scale_factor = 1.0;

    auto& t   = current_time();
    auto& u_n = current_solution();

    while (Dune::XT::Common::FloatCmp::gt(mixed_error, tol_)) {
      actual_dt *= time_step_scale_factor;

      for (size_t ii = 0; ii < num_stages_; ++ii) {
        // GRK4A and GRK4T use the same coefficients for the third and fourth stage
        if ((method == RosenbrockTimeStepperMethods::GRK4A || method == RosenbrockTimeStepperMethods::GRK4T) && ii == 3)
          u_intermediate_stages_[ii].vector() = u_intermediate_stages_[ii - 1].vector();
        else {
          u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
          u_tmp_.vector() = u_n.vector();
          for (size_t jj = 0; jj < ii; ++jj)
            u_tmp_.vector() += u_intermediate_stages_[jj].vector() * A_[ii][jj];
          op_.apply(u_tmp_, u_intermediate_stages_[ii], t + actual_dt * c_[ii]);
          u_intermediate_stages_[ii].vector() *= r_;
        }

        // TODO: If ii == 0, calculate (negative of) jacobian of L and partielle_Ableitung_von_L_nach_t
        //
        //
        //

        // if gamma is the same for all i, we only need to calculate the matrix in the first step
        if ((gamma_ii_equal_for_all_i_ && ii == 0) || !gamma_ii_equal_for_all_i_) {
          // create solver
          system_matrix_ = jacobian_ * -1.0 * r_;
          for (size_t row = 0; row < system_matrix_.rows(); ++row)
            system_matrix_.add_to_entry(row, row, 1.0 / (Gamma_[ii][ii] * actual_dt));
          solver_ = Dune::XT::Common::make_unique<SolverType>(system_matrix_);
        }
        for (size_t jj = 0; jj < ii; ++jj)
          u_intermediate_stages_[ii].vector().axpy(C_[ii][jj] / actual_dt, u_intermediate_stages_[jj].vector());
        u_intermediate_stages_[ii].vector() += partielle_Ableitung_von_L_nach_t_ * (r_ * actual_dt * d_[ii]);
        // solve, TODO: maybe need to do a copy before applying solver?
        solver_->apply(u_intermediate_stages_[ii].vector(), u_intermediate_stages_[ii].vector());
      }

      // compute error vector
      u_tmp_.vector() = u_intermediate_stages_[0].vector() * m_diff_[0];
      for (size_t ii = 1; ii < num_stages_; ++ii)
        u_tmp_.vector() += u_intermediate_stages_[ii].vector() * m_diff_[ii];

      // calculate u at timestep n+1
      for (size_t ii = 0; ii < num_stages_; ++ii)
        u_n.vector() += u_intermediate_stages_[ii].vector() * (m_1_[ii]);

      // scale error, use absolute error if norm is less than 0.01 and relative error else
      auto& diff_vector = u_tmp_.vector();
      for (size_t ii = 0; ii < diff_vector.size(); ++ii) {
        if (std::abs(u_n.vector()[ii]) > 0.01)
          diff_vector[ii] /= std::abs(u_n.vector()[ii]);
      }
      mixed_error = diff_vector.sup_norm();
      // scale dt to get the estimated optimal time step length, TODO: adapt formula
      time_step_scale_factor =
          std::min(std::max(0.9 * std::pow(tol_ / mixed_error, 1.0 / 4.0), scale_factor_min_), scale_factor_max_);

      if (mixed_error > tol_) { // go back from u at timestep n+1 to timestep n
        for (size_t ii = 0; ii < num_stages_; ++ii)
          u_n.vector() += u_intermediate_stages_[ii].vector() * (-1.0 * m_1_[ii]);
      }
    } // while (mixed_error > tol_)

    t += actual_dt;

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

private:
  const OperatorType& op_;
  const RangeFieldType r_;
  const RangeFieldType tol_;
  const TimeFieldType scale_factor_min_;
  const TimeFieldType scale_factor_max_;
  DiscreteFunctionType u_tmp_;
  const MatrixType A_;
  const VectorType m_1_;
  const VectorType m_2_;
  const VectorType c_;
  const VectorType m_diff_;
  VectorType d_;
  MatrixType C_;
  std::vector<DiscreteFunctionType> u_intermediate_stages_;
  typename SolverType::MatrixType system_matrix_, jacobian_;
  typename DiscreteFunctionType::VectorType partielle_Ableitung_von_L_nach_t_;
  const size_t num_stages_;
  std::unique_ptr<SolverType> solver_;
  const MatrixType Gamma_;
  bool gamma_ii_equal_for_all_i_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_TIMESTEPPER_ROSENBROCK_HH
