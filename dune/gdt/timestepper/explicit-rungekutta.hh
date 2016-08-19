// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/stuff/la/container.hh>

#include "interface.hh"


namespace Dune {
namespace GDT {


namespace internal {


// unspecialized
template <class RangeFieldType, class TimeFieldType, TimeStepperMethods method>
struct ButcherArrayProvider
{
  static_assert(AlwaysFalse<RangeFieldType>::value,
                "You cannot use ExplicitRungeKuttaTimeStepper with this value of TimeStepperMethods!");
};

// user-provided Butcher array
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, TimeStepperMethods::explicit_rungekutta_other>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicMatrix<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<TimeFieldType>();
  }
};

// Euler
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, TimeStepperMethods::explicit_euler>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0]");
  }
};

// Second order SSP
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, TimeStepperMethods::explicit_rungekutta_second_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0; 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5 0.5]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0 1]");
  }
};

// Third order SSP
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, TimeStepperMethods::explicit_rungekutta_third_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0 0; 1 0 0; 0.25 0.25 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " " + Dune::XT::Common::to_string(2.0 / 3.0, 15)
        + "]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0 1 0.5]");
  }
};

// Classic fourth order RK
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, TimeStepperMethods::explicit_rungekutta_classic_fourth_order>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 3.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 3.0, 15)
        + " "
        + Dune::XT::Common::to_string(1.0 / 6.0, 15)
        + "]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<TimeFieldType>>("[0 0.5 0.5 1]");
  }
};


} // namespace internal


/** \brief Time stepper using Runge Kutta methods
 *
 * Timestepper using explicit Runge Kutta methods to solve equations of the form u_t = r * L(u, t) where u is a
 * discrete function, L an operator acting on u and r a scalar factor (e.g. -1).
 * The specific Runge Kutta method can be chosen as the third template argument. If your desired Runge Kutta method is
 * not contained in ExplicitRungeKuttaMethods, choose ExplicitRungeKuttaMethods::other and supply a
 * DynamicMatrix< RangeFieldType > A and vectors (DynamicVector< RangeFieldType >) b and c in the constructor. Here, A,
 * b and c form the butcher tableau (see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods, A is
 * composed of the coefficients a_{ij}, b of b_j and c of c_j). The default is a forward euler method.
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values
 */
template <class OperatorImp, class DiscreteFunctionImp, class TimeFieldImp = double,
          TimeStepperMethods method = TimeStepperMethods::explicit_euler>
class ExplicitRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp> BaseType;
  typedef typename internal::ButcherArrayProvider<typename BaseType::RangeFieldType, TimeFieldImp, method>
      ButcherArrayProviderType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;

  typedef OperatorImp OperatorType;
  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::DynamicVector<RangeFieldType> VectorType;
  typedef typename Dune::DynamicVector<TimeFieldType> TimeVectorType;

  using BaseType::current_solution;
  using BaseType::current_time;

  /**
   * \brief Constructor for RungeKutta time stepper
   * \param op Operator L
   * \param initial_values Discrete function containing initial values for u at time t_0.
   * \param r Scalar factor (see above, default is 1)
   * \param t_0 Initial time (default is 0)
   * \param A Coefficient matrix (only provide if you use ExplicitRungeKuttaMethods::other)
   * \param b Coefficient vector (only provide if you use ExplicitRungeKuttaMethods::other)
   * \param c Coefficients for time steps (only provide if you use ExplicitRungeKuttaMethods::other)
   */
  ExplicitRungeKuttaTimeStepper(const OperatorType& op, const DiscreteFunctionType& initial_values,
                                const RangeFieldType r = 1.0, const double t_0 = 0.0,
                                const MatrixType& A     = ButcherArrayProviderType::A(),
                                const VectorType& b     = ButcherArrayProviderType::b(),
                                const TimeVectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , op_(op)
    , r_(r)
    , u_tmp_(BaseType::current_solution())
    , A_(A)
    , b_(b)
    , c_(c)
    , num_stages_(A_.rows())
  {
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_.size() == A_.rows());
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
      u_intermediate_stages_.emplace_back(current_solution());
    }
  } // constructor

  /**
   * \brief Constructor ignoring the tol argument for compatibility with AdaptiveRungeKuttaTimeStepper
   */
  ExplicitRungeKuttaTimeStepper(const OperatorType& op, const DiscreteFunctionType& initial_values,
                                const RangeFieldType r, const double t_0, const RangeFieldType /*tol*/)
    : ExplicitRungeKuttaTimeStepper(op, initial_values, r, t_0)
  {
  }

  virtual TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt) override final
  {
    const TimeFieldType actual_dt = std::min(dt, max_dt);
    auto& t                       = current_time();
    auto& u_n                     = current_solution();
    // calculate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
      u_tmp_.vector() = u_n.vector();
      for (size_t jj = 0; jj < ii; ++jj)
        u_tmp_.vector() += u_intermediate_stages_[jj].vector() * (actual_dt * r_ * (A_[ii][jj]));
      op_.apply(u_tmp_, u_intermediate_stages_[ii], t + actual_dt * c_[ii]);
    }

    // calculate value of u at next time step
    for (size_t ii = 0; ii < num_stages_; ++ii)
      u_n.vector() += u_intermediate_stages_[ii].vector() * (r_ * actual_dt * b_[ii]);

    // augment time
    t += actual_dt;

    return dt;
  } // ... step(...)

  const std::pair<bool, TimeFieldType>
  find_suitable_dt(const TimeFieldType initial_dt, const TimeFieldType dt_refinement_factor = 2,
                   const RangeFieldType treshold = 0.9 * std::numeric_limits<RangeFieldType>::max(),
                   const size_t max_steps_per_dt = 20, const size_t max_refinements = 20)
  {
    auto& t   = current_time();
    auto& u_n = current_solution();
    assert(treshold > 0);
    // save current state
    DiscreteFunctionType initial_u_n = u_n;
    TimeFieldType initial_t          = t;
    // start with initial dt
    TimeFieldType current_dt = initial_dt;
    size_t num_refinements   = 0;
    while (num_refinements < max_refinements) {
      std::cout << "Trying time step length dt = " << current_dt << "... " << std::flush;
      bool unlikely_value_occured = false;
      size_t num_steps            = 0;
      // do max_steps_per_dt time steps...
      while (!unlikely_value_occured) {
        step(current_dt);
        ++num_steps;
        // ... unless there is a value above threshold
        for (size_t kk = 0; kk < u_n.vector().size(); ++kk) {
          if (std::abs(u_n.vector()[kk]) > treshold) {
            unlikely_value_occured = true;
            std::cout << "failed" << std::endl;
            break;
          }
        }
        // if we are able to do max_steps_per_dt time steps with this dt, we accept this dt
        if (num_steps == max_steps_per_dt) {
          std::cout << "looks fine" << std::endl;
          u_n.vector() = initial_u_n.vector();
          t            = initial_t;
          return std::make_pair(bool(true), current_dt);
        }
      }
      // if there was a value above threshold start over with smaller dt
      u_n.vector() = initial_u_n.vector();
      t            = initial_t;
      current_dt /= dt_refinement_factor;
      ++num_refinements;
    }
    return std::make_pair(bool(false), current_dt);
  }

private:
  const OperatorType& op_;
  const RangeFieldType r_;
  DiscreteFunctionType u_tmp_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector<DiscreteFunctionType> u_intermediate_stages_;
  const size_t num_stages_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH
