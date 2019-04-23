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

#ifndef DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH

#include <utility>

#include "enums.hh"
#include "interface.hh"


namespace Dune {
namespace GDT {


namespace internal {


// unspecialized
template <class RangeFieldType, TimeStepperMethods method>
struct ButcherArrayProvider
{
  static_assert(AlwaysFalse<RangeFieldType>::value,
                "You cannot use ExplicitRungeKuttaTimeStepper with this value of TimeStepperMethods!");
};

// user-provided Butcher array
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_other>
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

  static Dune::DynamicVector<RangeFieldType> c()
  {
    DUNE_THROW(Dune::NotImplemented,
               "You have to provide a Butcher array in ExplicitRungeKuttaTimeStepper's constructor for this method!");
    return Dune::DynamicVector<RangeFieldType>();
  }
};

// Euler
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_euler>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[1]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0]");
  }
};

// Second order SSP
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_second_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0; 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0.5 0.5]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0 1]");
  }
};

// Third order SSP
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_third_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>("[0 0 0; 1 0 0; 0.25 0.25 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " "
        + Dune::XT::Common::to_string(2.0 / 3.0, 15) + "]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0 1 0.5]");
  }
};

// Classic fourth order RK
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_classic_fourth_order>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return Dune::XT::Common::from_string<Dune::DynamicMatrix<RangeFieldType>>(
        "[0 0 0 0; 0.5 0 0 0; 0 0.5 0 0; 0 0 1 0]");
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>(
        "[" + Dune::XT::Common::to_string(1.0 / 6.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 3.0, 15) + " "
        + Dune::XT::Common::to_string(1.0 / 3.0, 15) + " " + Dune::XT::Common::to_string(1.0 / 6.0, 15) + "]");
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return Dune::XT::Common::from_string<Dune::DynamicVector<RangeFieldType>>("[0 0.5 0.5 1]");
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
 * Notation: For an s-stage method,
 * \mathbf{u}^{n+1} = \mathbf{u}^n + dt \sum_{i=0}^{s-1} b_i \mathbf{k}_i
 * \mathbf{k}_i = L(\mathbf{u}_i, t^n + dt c_i)
 * \mathbf{u}_i = \mathbf{u}^n + dt \sum_{j=0}^{i-1} a_{ij} \mathbf{k}_j,
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values
 */
template <class OperatorImp, class DiscreteFunctionImp, TimeStepperMethods method = TimeStepperMethods::explicit_euler>
class ExplicitRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  using BaseType = TimeStepperInterface<DiscreteFunctionImp>;
  using ButcherArrayProviderType = typename internal::ButcherArrayProvider<typename BaseType::RangeFieldType, method>;

public:
  using typename BaseType::DataHandleType;
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DiscreteSolutionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;

  using OperatorType = OperatorImp;
  using MatrixType = Dune::DynamicMatrix<RangeFieldType>;
  using VectorType = Dune::DynamicVector<RangeFieldType>;

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
  ExplicitRungeKuttaTimeStepper(const OperatorType& op,
                                DiscreteFunctionType& initial_values,
                                const RangeFieldType r = 1.0,
                                const double t_0 = 0.0,
                                const MatrixType& A = ButcherArrayProviderType::A(),
                                const VectorType& b = ButcherArrayProviderType::b(),
                                const VectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , op_(&op)
    , r_(r)
    , u_i_(BaseType::current_solution())
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
    // store as many discrete functions as needed for the stages k
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(current_solution());
    }
  } // constructor

  /**
   * \brief Constructor ignoring the tol argument for compatibility with AdaptiveRungeKuttaTimeStepper
   */
  ExplicitRungeKuttaTimeStepper(const OperatorType& op,
                                const DiscreteFunctionType& initial_values,
                                const RangeFieldType r,
                                const double t_0,
                                const RangeFieldType /*tol*/)
    : ExplicitRungeKuttaTimeStepper(op, initial_values, r, t_0)
  {}

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    // calculate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_i_.dofs().vector() = u_n.dofs().vector();
      for (size_t jj = 0; jj < ii; ++jj)
        u_i_.dofs().vector() += stages_k_[jj].dofs().vector() * (actual_dt * r_ * (A_[ii][jj]));
      // TODO: provide actual_dt to op_. This leads to spurious oscillations in the Lax-Friedrichs flux
      // because actual_dt/dx may become very small.
      op_->apply(u_i_.dofs().vector(),
                 stages_k_[ii].dofs().vector(),
                 XT::Common::Parameter({{"t", {t + actual_dt * c_[ii]}}, {"dt", {dt}}}));
      DataHandleType stages_k_ii_handle(stages_k_[ii]);
      stages_k_[ii].space().grid_view().template communicate<DataHandleType>(
          stages_k_ii_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    }

    // calculate value of u at next time step
    for (size_t ii = 0; ii < num_stages_; ++ii)
      u_n.dofs().vector() += stages_k_[ii].dofs().vector() * (r_ * actual_dt * b_[ii]);

    // augment time
    t += actual_dt;

    return dt;
  } // ... step(...)


  void set_operator(const OperatorType& op)
  {
    op_ = &op;
  }

  const std::pair<bool, RangeFieldType>
  find_suitable_dt(const RangeFieldType initial_dt,
                   const RangeFieldType dt_refinement_factor = 2,
                   const RangeFieldType treshold = 0.9 * std::numeric_limits<RangeFieldType>::max(),
                   const size_t max_steps_per_dt = 20,
                   const size_t max_refinements = 20)
  {
    auto& t = current_time();
    auto& u_n = current_solution();
    assert(treshold > 0);
    // save current state
    DiscreteFunctionType initial_u_n = u_n;
    RangeFieldType initial_t = t;
    // start with initial dt
    RangeFieldType current_dt = initial_dt;
    size_t num_refinements = 0;
    while (num_refinements < max_refinements) {
      std::cout << "Trying time step length dt = " << current_dt << "... " << std::flush;
      bool unlikely_value_occured = false;
      size_t num_steps = 0;
      // do max_steps_per_dt time steps...
      while (!unlikely_value_occured) {
        step(current_dt);
        ++num_steps;
        // ... unless there is a value above threshold
        for (size_t kk = 0; kk < u_n.dofs().vector().size(); ++kk) {
          if (std::abs(u_n.dofs().vector()[kk]) > treshold) {
            unlikely_value_occured = true;
            std::cout << "failed" << std::endl;
            break;
          }
        }
        // if we are able to do max_steps_per_dt time steps with this dt, we accept this dt
        if (num_steps == max_steps_per_dt) {
          std::cout << "looks fine" << std::endl;
          u_n.dofs().vector() = initial_u_n.dofs().vector();
          t = initial_t;
          return std::make_pair(bool(true), current_dt);
        }
      }
      // if there was a value above threshold start over with smaller dt
      u_n.dofs().vector() = initial_u_n.dofs().vector();
      t = initial_t;
      current_dt /= dt_refinement_factor;
      ++num_refinements;
    }
    return std::make_pair(bool(false), current_dt);
  }

private:
  const OperatorType* op_;
  const RangeFieldType r_;
  DiscreteFunctionType u_i_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector<DiscreteFunctionType> stages_k_;
  const size_t num_stages_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH
