// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_TIMESTEPPER_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_RUNGEKUTTA_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/la/container.hh>

#include "interface.hh"


namespace Dune {
namespace GDT {
namespace TimeStepper {


enum class RungeKuttaMethods
{
  euler,
  second_order_ssp,
  third_order_ssp,
  classic_fourth_order,
  other
};


namespace internal {


// unspecialized
template <class RangeFieldType, class TimeFieldType, RungeKuttaMethods method = RungeKuttaMethods::other>
struct ButcherArrayProvider
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttas constructor for this method!");
    return Dune::DynamicMatrix<RangeFieldType>();
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttas constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttas constructor for this method!");
    return Dune::DynamicVector<TimeFieldType>();
  }
};


// Euler
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, RungeKuttaMethods::euler>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return DSC::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return DSC::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return DSC::from_string<Dune::DynamicVector<TimeFieldType>>("[0]");
  }
};

// Second order SSP
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, RungeKuttaMethods::second_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return DSC::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0; 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return DSC::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5 0.5]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return DSC::from_string<Dune::DynamicVector<TimeFieldType>>("[0 1]");
  }
};

// Third order SSP
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, RungeKuttaMethods::third_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return DSC::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0 0; 1 0 0; 0.25 0.25 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return DSC::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + DSC::to_string(1.0 / 6.0, 15) + " " + DSC::to_string(1.0 / 6.0, 15) + " " + DSC::to_string(2.0 / 3.0, 15)
        + "]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return DSC::from_string<Dune::DynamicVector<TimeFieldType>>("[0 1 0.5]");
  }
};

// Classic fourth order RK
template <class RangeFieldType, class TimeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeFieldType, RungeKuttaMethods::classic_fourth_order>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return DSC::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return DSC::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + DSC::to_string(1.0 / 6.0, 15) + " " + DSC::to_string(1.0 / 3.0, 15) + " " + DSC::to_string(1.0 / 3.0, 15)
        + " "
        + DSC::to_string(1.0 / 6.0, 15)
        + "]");
  }

  static Dune::DynamicVector<TimeFieldType> c()
  {
    return DSC::from_string<Dune::DynamicVector<TimeFieldType>>("[0 0.5 0.5 1]");
  }
};


} // namespace internal


/** \brief Time stepper using Runge Kutta methods
 *
 * Timestepper using explicit Runge Kutta methods to solve equations of the form u_t = \alpha L(u) where u is a
 * discrete function, L an operator acting on u and \alpha a scalar factor (e.g. -1).
 * The specific Runge Kutta method can be chosen as the third template argument. If your desired Runge Kutta method is
 * not contained in Dune::GDT::TimeStepper::RungeKuttaMethods, choose RungeKuttaMethods::other and supply a
 * DynamicMatrix< RangeFieldType > A and vectors (DynamicVector< RangeFieldType >) b and c in the constructor. Here, A,
 * b and c form the butcher tableau (see https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods, A is
 * composed of the coefficients a_{ij}, b of b_j and c of c_j). The default is a forward euler method.
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values
 */
template <class OperatorImp, class DiscreteFunctionImp, class TimeFieldImp,
          RungeKuttaMethods method = RungeKuttaMethods::euler>
class ExplicitRungeKutta : public TimeStepperInterface<DiscreteFunctionImp, TimeFieldImp>
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

  /**
   * \brief Constructor for RungeKutta time stepper
   *
   * \param flux operator L
   * \param initial_values Discrete function containing initial values for u
   * \param A A (see above)
   * \param b b (see above)
   * \param c c (completely ignored, see above)
   */
  ExplicitRungeKutta(const OperatorType& op, const DiscreteFunctionType& initial_values,
                     const RangeFieldType alpha = 1.0, const double t_0 = 0.0,
                     const MatrixType A     = ButcherArrayProviderType::A(),
                     const VectorType b     = ButcherArrayProviderType::b(),
                     const TimeVectorType c = ButcherArrayProviderType::c())
    : op_(op)
    , initial_values_(initial_values)
    , alpha_(alpha)
    , u_n_(initial_values_)
    , u_tmp_(u_n_)
    , t_(t_0)
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
        assert(DSC::FloatCmp::eq(A_[ii][jj], 0.0)
               && "A has to be a lower triangular matrix with 0 on the main diagonal");
      }
    }
#endif // NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_intermediate_stages_.emplace_back(u_n_);
    }
  } // constructor

  virtual TimeFieldType current_time() const override
  {
    return t_;
  }

  virtual const DiscreteFunctionType& solution_at_current_time() const override
  {
    return u_n_;
  }

  virtual void set_solution_at_current_time(const DiscreteFunctionType& discrete_func) override final
  {
    u_n_.vector() = discrete_func.vector();
  }

  TimeFieldType step(const TimeFieldType dt) override final
  {
    // calculate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
      u_tmp_.vector() = u_n_.vector();
      for (size_t jj = 0; jj < ii; ++jj)
        u_tmp_.vector() += u_intermediate_stages_[jj].vector() * (dt * alpha_ * (A_[ii][jj]));
      op_.apply(u_tmp_, u_intermediate_stages_[ii], t_ + dt * c_[ii]);
    }

    // calculate value of u at next time step
    for (size_t ii = 0; ii < num_stages_; ++ii)
      u_n_.vector() += u_intermediate_stages_[ii].vector() * (alpha_ * dt * b_[ii]);

    // augment time
    t_ += dt;

    return dt;
  } // ... step(...)

  const std::pair<bool, TimeFieldType>
  find_suitable_dt(const TimeFieldType initial_dt, const TimeFieldType dt_refinement_factor = 2,
                   const RangeFieldType treshold = 0.9 * std::numeric_limits<RangeFieldType>::max(),
                   const size_t max_steps_per_dt = 20, const size_t max_refinements = 20)
  {
    assert(treshold > 0);
    // save current state
    DiscreteFunctionType initial_u_n = u_n_;
    TimeFieldType initial_t          = t_;
    // start with initial dt
    TimeFieldType current_dt = initial_dt;
    size_t num_refinements = 0;
    while (num_refinements < max_refinements) {
      std::cout << "Trying time step length dt = " << current_dt << "... " << std::flush;
      bool unlikely_value_occured = false;
      size_t num_steps            = 0;
      // do max_steps_per_dt time steps...
      while (!unlikely_value_occured) {
        step(current_dt);
        ++num_steps;
        // ... unless there is a value above threshold
        for (size_t kk = 0; kk < u_n_.vector().size(); ++kk) {
          if (std::abs(u_n_.vector()[kk]) > treshold) {
            unlikely_value_occured = true;
            std::cout << "failed" << std::endl;
            break;
          }
        }
        // if we are able to do max_steps_per_dt time steps with this dt, we accept this dt
        if (num_steps == max_steps_per_dt) {
          std::cout << "looks fine" << std::endl;
          u_n_.vector() = initial_u_n.vector();
          t_ = initial_t;
          return std::make_pair(bool(true), current_dt);
        }
      }
      // if there was a value above threshold start over with smaller dt
      u_n_.vector() = initial_u_n.vector();
      t_ = initial_t;
      current_dt /= dt_refinement_factor;
      ++num_refinements;
    }
    return std::make_pair(bool(false), current_dt);
  }

private:
  const OperatorType& op_;
  const DiscreteFunctionType& initial_values_;
  const RangeFieldType alpha_;
  DiscreteFunctionType u_n_;
  DiscreteFunctionType u_tmp_;
  double t_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector<DiscreteFunctionType> u_intermediate_stages_;
  const size_t num_stages_;
};


template <class FirstStepperImp, class SecondStepperImp>
class FractionalStepStepper : public TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType,
                                                          typename FirstStepperImp::TimeFieldType>
{
  typedef TimeStepperInterface<typename FirstStepperImp::DiscreteFunctionType, typename FirstStepperImp::TimeFieldType>
      BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::TimeFieldType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::SolutionType;

  typedef FirstStepperImp FirstStepperType;
  typedef SecondStepperImp SecondStepperType;

  FractionalStepStepper(FirstStepperType& first_stepper, SecondStepperType& second_stepper)
    : first_stepper_(first_stepper)
    , second_stepper_(second_stepper)
  {
  } // constructor

  virtual TimeFieldType current_time() const override
  {
    return first_stepper_.current_time();
  }

  virtual const DiscreteFunctionType& solution_at_current_time() const override
  {
    return first_stepper_.solution_at_current_time();
  }

  virtual void set_solution_at_current_time(const DiscreteFunctionType& discrete_func) override final
  {
    first_stepper_.set_solution_at_current_time(discrete_func);
  }

  TimeFieldType step(const TimeFieldType dt) override final
  {
    const auto t = current_time();
    first_stepper_.solve(t + dt, dt, -1, false);
    second_stepper_.set_solution_at_current_time(first_stepper_.solution_at_current_time());
    second_stepper_.solve(t + dt, dt, -1, false);
    first_stepper_.set_solution_at_current_time(second_stepper_.solution_at_current_time());
    return dt;
  } // ... step(...)

private:
  FirstStepperType& first_stepper_;
  SecondStepperType& second_stepper_;
};


} // namespace TimeStepper
} // namespace Stuff
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_RUNGEKUTTA_HH
