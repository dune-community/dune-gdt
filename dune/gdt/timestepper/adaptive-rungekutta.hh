// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TIMESTEPPER_ADAPTIVE_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_ADAPTIVE_RUNGEKUTTA_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include "interface.hh"


namespace Dune {
namespace GDT {


namespace internal {


// unspecialized
template <class RangeFieldType, TimeStepperMethods method>
struct AdaptiveButcherArrayProvider
{
  static_assert(AlwaysFalse<RangeFieldType>::value,
                "You cannot use AdaptiveRungeKuttaTimeStepper with this value of TimeStepperMethods!");
};

// user-provided Butcher array
template <class RangeFieldType>
struct AdaptiveButcherArrayProvider<RangeFieldType, TimeStepperMethods::adaptive_rungekutta_other>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in AdaptiveRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicMatrix<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in AdaptiveRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in AdaptiveRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in AdaptiveRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }
};

// Bogacki-Shampine (adaptive RK23)
template <class RangeFieldType>
class AdaptiveButcherArrayProvider<RangeFieldType, TimeStepperMethods::bogacki_shampine>
{
public:
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(
        "[0 0 0 0; 0.5 0 0 0; 0 0.75 0 0; " + Dune::XT::Common::to_string(2.0 / 9.0, 15) + " "
        + Dune::XT::Common::to_string(1.0 / 3.0, 15)
        + " "
        + Dune::XT::Common::to_string(4.0 / 9.0, 15)
        + " 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(2.0 / 9.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 3.0, 15) + " "
        + Dune::XT::Common::to_string(4.0 / 9.0, 15)
        + " 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(7.0 / 24.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 4.0, 15) + " "
        + Dune::XT::Common::to_string(1.0 / 3.0, 15)
        + " "
        + Dune::XT::Common::to_string(1.0 / 8.0, 15)
        + " 0]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5 0.75 1 0]");
  }
};

// Dormand-Prince (adaptive RK45)
template <class RangeFieldType>
class AdaptiveButcherArrayProvider<RangeFieldType, TimeStepperMethods::dormand_prince>
{
public:
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(
        std::string("[0 0 0 0 0 0 0;") + " 0.2 0 0 0 0 0 0;" + " 0.075 0.225 0 0 0 0 0;" + " "
        + Dune::XT::Common::to_string(44.0 / 45.0, 15)
        + " "
        + Dune::XT::Common::to_string(-56.0 / 15.0, 15)
        + " "
        + Dune::XT::Common::to_string(32.0 / 9.0, 15)
        + " 0 0 0 0;"
        + " "
        + Dune::XT::Common::to_string(19372.0 / 6561.0, 15)
        + " "
        + Dune::XT::Common::to_string(-25360.0 / 2187.0, 15)
        + " "
        + Dune::XT::Common::to_string(64448.0 / 6561.0, 15)
        + " "
        + Dune::XT::Common::to_string(-212.0 / 729.0, 15)
        + " 0 0 0;"
        + " "
        + Dune::XT::Common::to_string(9017.0 / 3168.0, 15)
        + " "
        + Dune::XT::Common::to_string(-355.0 / 33.0, 15)
        + " "
        + Dune::XT::Common::to_string(46732.0 / 5247.0, 15)
        + " "
        + Dune::XT::Common::to_string(49.0 / 176.0, 15)
        + " "
        + Dune::XT::Common::to_string(-5103.0 / 18656.0, 15)
        + " 0 0;"
        + " "
        + Dune::XT::Common::to_string(35.0 / 384.0, 15)
        + " 0 "
        + Dune::XT::Common::to_string(500.0 / 1113.0, 15)
        + " "
        + Dune::XT::Common::to_string(125.0 / 192.0, 15)
        + " "
        + Dune::XT::Common::to_string(-2187.0 / 6784.0, 15)
        + " "
        + Dune::XT::Common::to_string(11.0 / 84.0, 15)
        + " 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_1()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(35.0 / 384.0, 15) + " 0 " + Dune::XT::Common::to_string(500.0 / 1113.0, 15)
        + " "
        + Dune::XT::Common::to_string(125.0 / 192.0, 15)
        + " "
        + Dune::XT::Common::to_string(-2187.0 / 6784.0, 15)
        + " "
        + Dune::XT::Common::to_string(11.0 / 84.0, 15)
        + " 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b_2()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(5179.0 / 57600.0, 15) + " 0 "
        + Dune::XT::Common::to_string(7571.0 / 16695.0, 15)
        + " "
        + Dune::XT::Common::to_string(393.0 / 640.0, 15)
        + " "
        + Dune::XT::Common::to_string(-92097.0 / 339200.0, 15)
        + " "
        + Dune::XT::Common::to_string(187.0 / 2100.0, 15)
        + " "
        + Dune::XT::Common::to_string(1.0 / 40.0, 15)
        + "]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[0 0.2 0.3 0.8 " + Dune::XT::Common::to_string(8.0 / 9.0, 15) + " 1 1]");
  }
}; // Dormand-Prince (RK45)


} // namespace internal


/** \brief Time stepper using adaptive Runge Kutta methods
 *
 * Timestepper using adaptive Runge Kutta methods to solve equations of the form u_t = r * L(u, t) where u is a
 * discrete function, L an operator acting on u and \alpha a scalar factor (e.g. -1).
 * The specific Runge Kutta method can be chosen as the third template argument. If your desired Runge Kutta method is
 * not contained in AdaptiveRungeKuttaMethods, choose AdaptiveRungeKuttaMethods::other and
 * supply a DynamicMatrix< RangeFieldType > A and vectors b_1, b_2 (DynamicVector< RangeFieldType >) and c
 * (DynamicVector< RangeFieldType >) in the constructor. Here, A, b_1, b_2 and c form the butcher tableau (see
 * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Embedded_methods, A is composed of the coefficients
 * a_{ij}, b_1 of b_j, b_2 of b_j^* and c of c_j). The default is the Dormand-Prince RK45 method.
 * In each time step, the error is estimated using the difference between the two solutions obtained using either b_1 or
 * b_2. If the estimated error is higher than a specified tolerance tol, the calculation is repeated with a smaller
 * time step. The tolerance tol and the error estimate are also used to estimate the optimal time step length for the
 * next time step via dt_new = dt_old*min(max(0.9*(tol/error)^(1/5), scale_factor_min), scale_factor_max_);
 *
 * Notation: For an s-stage method,
 * \mathbf{u}^{n+1} = \mathbf{u}^n + dt \sum_{i=0}^{s-1} b_i \mathbf{k}_i
 * \mathbf{k}_i = L(\mathbf{u}_i, t^n + dt c_i)
 * \mathbf{u}_i = \mathbf{u}^n + dt \sum_{j=0}^{i-1} a_{ij} \mathbf{k}_j,
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values and solution at a fixed time
 * \tparam method Adaptive Runge-Kutta method that is used (default is AdaptiveRungeKuttaMethods::dormand_prince)
 */
template <class OperatorImp, class DiscreteFunctionImp, TimeStepperMethods method = TimeStepperMethods::dormand_prince>
class AdaptiveRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp> BaseType;
  typedef typename internal::AdaptiveButcherArrayProvider<typename BaseType::RangeFieldType, method>
      ButcherArrayProviderType;

public:
  typedef OperatorImp OperatorType;
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::DynamicVector<RangeFieldType> VectorType;
  typedef typename std::vector<std::pair<RangeFieldType, DiscreteFunctionType>> SolutionType;

  /**
   * \brief Constructor for AdaptiveRungeKuttaTimeStepper time stepper
   * \param op Operator L
   * \param initial_values Discrete function containing initial values for u at time t_0.
   * \param r Scalar factor (see above, default is 1)
   * \param t_0 Initial time (default is 0)
   * \param tol Error tolerance for the adaptive scheme (default is 1e-4)
   * \param scale_factor_min Minimum allowed factor for time step scaling (default is 0.2).
   * \param scale_factor_max Maximum allowed factor for time step scaling (default is 5).
   * \param A Coefficient matrix (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param b_1 First set of coefficients (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param b_2 Second set of coefficients (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param c Coefficients for time steps (only provide if you use AdaptiveRungeKuttaMethods::other)
   */
  AdaptiveRungeKuttaTimeStepper(const OperatorType& op,
                                DiscreteFunctionType& initial_values,
                                const RangeFieldType r = 1.0,
                                const double t_0 = 0.0,
                                const RangeFieldType tol = 1e-4,
                                const RangeFieldType scale_factor_min = 0.2,
                                const RangeFieldType scale_factor_max = 5,
                                const MatrixType& A = ButcherArrayProviderType::A(),
                                const VectorType& b_1 = ButcherArrayProviderType::b_1(),
                                const VectorType& b_2 = ButcherArrayProviderType::b_2(),
                                const VectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , op_(op)
    , r_(r)
    , tol_(tol)
    , scale_factor_min_(scale_factor_min)
    , scale_factor_max_(scale_factor_max)
    , u_tmp_(BaseType::current_solution())
    , A_(A)
    , b_1_(b_1)
    , b_2_(b_2)
    , c_(c)
    , b_diff_(b_2_ - b_1_)
    , num_stages_(A_.rows())
  {
    assert(Dune::XT::Common::FloatCmp::gt(tol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::le(scale_factor_min_, 1.0));
    assert(Dune::XT::Common::FloatCmp::ge(scale_factor_max_, 1.0));
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_1_.size() == A_.rows());
    assert(b_2_.size() == A_.rows());
    assert(c_.size() == A_.rows());
#ifndef NDEBUG
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        assert(Dune::XT::Common::FloatCmp::eq(A_[ii][jj], 0.0)
               && "A has to be a lower triangular matrix with 0 on the main diagonal");
      }
    }
#endif // NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(current_solution());
    }
  } // constructor AdaptiveRungeKuttaTimeStepper

  using BaseType::current_solution;
  using BaseType::current_time;

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt)
  {
    RangeFieldType actual_dt = std::min(dt, max_dt);
    RangeFieldType mixed_error = std::numeric_limits<RangeFieldType>::max();
    RangeFieldType time_step_scale_factor = 1.0;

    auto& t = current_time();
    auto& u_n = current_solution();

    while (Dune::XT::Common::FloatCmp::gt(mixed_error, tol_)) {
      actual_dt *= time_step_scale_factor;
      size_t first_stage_to_compute = 0;
      if (last_stage_of_previous_step_) {
        stages_k_[0].vector() = last_stage_of_previous_step_->vector();
        first_stage_to_compute = 1;
      }

      for (size_t ii = first_stage_to_compute; ii < num_stages_; ++ii) {
        stages_k_[ii].vector() *= RangeFieldType(0);
        u_tmp_.vector() = u_n.vector();
        for (size_t jj = 0; jj < ii; ++jj)
          u_tmp_.vector() += stages_k_[jj].vector() * (actual_dt * r_ * (A_[ii][jj]));
        op_.apply(u_tmp_, stages_k_[ii], t + actual_dt * c_[ii]);
      }

      // compute error vector
      u_tmp_.vector() = stages_k_[0].vector() * b_diff_[0];
      for (size_t ii = 1; ii < num_stages_; ++ii)
        u_tmp_.vector() += stages_k_[ii].vector() * b_diff_[ii];
      u_tmp_.vector() *= actual_dt * r_;

      // calculate u at timestep n+1
      for (size_t ii = 0; ii < num_stages_; ++ii)
        u_n.vector() += stages_k_[ii].vector() * (actual_dt * r_ * b_1_[ii]);

      // scale error, use absolute error if norm is less than 0.01 and relative error else
      auto& diff_vector = u_tmp_.vector();
      for (size_t ii = 0; ii < diff_vector.size(); ++ii) {
        if (std::abs(u_n.vector()[ii]) > 0.01)
          diff_vector[ii] /= std::abs(u_n.vector()[ii]);
      }
      mixed_error = diff_vector.sup_norm();
      // scale dt to get the estimated optimal time step length, TODO: adapt formula
      time_step_scale_factor =
          std::min(std::max(0.9 * std::pow(tol_ / mixed_error, 1.0 / 5.0), scale_factor_min_), scale_factor_max_);

      if (mixed_error > tol_) { // go back from u at timestep n+1 to timestep n
        for (size_t ii = 0; ii < num_stages_; ++ii)
          u_n.vector() += stages_k_[ii].vector() * (-1.0 * r_ * actual_dt * b_1_[ii]);
      }
    } // while (mixed_error > tol_)

    if (!last_stage_of_previous_step_)
      last_stage_of_previous_step_ = Dune::XT::Common::make_unique<DiscreteFunctionType>(u_n);
    last_stage_of_previous_step_->vector() = stages_k_[num_stages_ - 1].vector();

    t += actual_dt;

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

private:
  const OperatorType& op_;
  const RangeFieldType r_;
  const RangeFieldType tol_;
  const RangeFieldType scale_factor_min_;
  const RangeFieldType scale_factor_max_;
  DiscreteFunctionType u_tmp_;
  const MatrixType A_;
  const VectorType b_1_;
  const VectorType b_2_;
  const VectorType c_;
  const VectorType b_diff_;
  std::vector<DiscreteFunctionType> stages_k_;
  const size_t num_stages_;
  std::unique_ptr<DiscreteFunctionType> last_stage_of_previous_step_;
}; // class AdaptiveRungeKuttaTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_ADAPTIVE_RUNGEKUTTA_HH
