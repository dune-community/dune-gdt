// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_INTERFACE_HH
#define DUNE_GDT_TIMESTEPPER_INTERFACE_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/stuff/la/container.hh>


namespace Dune {
namespace GDT {


enum class TimeStepperMethods
{
  bogacki_shampine,
  dormand_prince,
  adaptive_rungekutta_other,
  explicit_euler,
  explicit_rungekutta_second_order_ssp,
  explicit_rungekutta_third_order_ssp,
  explicit_rungekutta_classic_fourth_order,
  explicit_rungekutta_other
};


namespace internal {


struct FloatCmpLt
{
  template <class A, class B>
  bool operator()(const A& a, const B& b) const
  {
    return Dune::XT::Common::FloatCmp::lt(a, b);
  }
};


} // namespace internal


template <class DiscreteFunctionImp, class TimeFieldImp>
class TimeStepperInterface
    : Dune::XT::Common::StorageProvider<DiscreteFunctionImp>,
      Dune::XT::Common::StorageProvider<std::map<TimeFieldImp, DiscreteFunctionImp, typename internal::FloatCmpLt>>
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef TimeFieldImp TimeFieldType;

  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename std::map<TimeFieldType, DiscreteFunctionType, typename internal::FloatCmpLt> SolutionType;
  typedef typename SolutionType::value_type TimeAndDiscreteFunctionPairType;

private:
  typedef typename Dune::XT::Common::StorageProvider<DiscreteFunctionImp> CurrentSolutionStorageProviderType;
  typedef typename Dune::XT::Common::StorageProvider<std::map<TimeFieldImp, DiscreteFunctionImp, typename internal::FloatCmpLt>>
      SolutionStorageProviderType;

protected:
  TimeStepperInterface(const TimeFieldType t_0, const DiscreteFunctionType& initial_values)
    : CurrentSolutionStorageProviderType(initial_values)
    , SolutionStorageProviderType(new SolutionType())
    , t_(t_0)
    , u_n_(&CurrentSolutionStorageProviderType::access())
    , solution_(&SolutionStorageProviderType::access())
  {
  }

public:
  virtual ~TimeStepperInterface() = default;

  /**
   * \brief Perform one time step and return the estimated optimal time step length for the next time step.
   * \param dt Time step length.
   * \param max_dt Maximal allowed time step length in this time step.
   * \return Estimated optimal time step length for the next time step.
   * \note For adaptive methods or if max_dt < dt, dt is just the initial time step length, the actual time step taken
   * may be shorter. To get the actual time step length you need to compare current_time() before and after calling
   * step().
   * \note This method has to increase the current time by the actual time step length taken, otherwise
   * solve will never finish execution (if current_time is not increased) or give wrong results (if current time is
   * increased by a wrong dt).
   */
  virtual TimeFieldType step(const TimeFieldType dt, const TimeFieldType max_dt) = 0;

  const TimeFieldType& current_time() const
  {
    return t_;
  }

  TimeFieldType& current_time()
  {
    return t_;
  }

  const DiscreteFunctionType& current_solution() const
  {
    return *u_n_;
  }

  DiscreteFunctionType& current_solution()
  {
    return *u_n_;
  }

  void set_current_solution_pointer(DiscreteFunctionType& discrete_function)
  {
    u_n_ = &discrete_function;
  }

  const SolutionType& solution() const
  {
    return *solution_;
  }

  SolutionType& solution()
  {
    return *solution_;
  }

  void set_solution_pointer(SolutionType& solution_ref)
  {
    solution_ = &solution_ref;
  }

  /**
   * \brief solve Applies time stepping scheme up to time t_end
   * \param initial_dt Initial time step length. Non-adaptive schemes will use this dt for all time steps. Adaptive
   * schemes will use initial_dt as a first guess for the time step length.
   * \param num_save_steps Number of time points that will be stored in solution if save_solution is true and/or
   * visualized if visualize_solution is true. Default is size_t(-1), which will store all time steps.
   * \param save_solution If true, the discrete solution at num_save_steps equidistant time points will be stored in
   * solution.
   * \param output_progress If true, the current time and current time step length will be written to std::cout at
   * num_save_steps equidistant time points
   * \param visualize If true, solution will be written to .vtp/.vtu files at num_save_steps equidistant time points
   * \param filename_prefix Prefix of .vtp/.vtu files for visualization
   * \param sol Optional object that solution will be appended to instead of storing it in the time stepper.
   * \return estimated optimal time step length for next step
   * \note If num_save_steps is specified (i.e. if it is not size_t(-1), the solution will be stored/visualized at
   * exactly num_save_steps + 1 equidistant time points (including the initial time and t_end), even if the time step
   * length has to be reduced to hit these time points.
   */
  virtual TimeFieldType solve(const TimeFieldType t_end, const TimeFieldType initial_dt, const size_t num_save_steps,
                              const bool save_solution, const bool output_progress, const bool visualize,
                              const std::string filename_prefix, SolutionType& sol)
  {
    TimeFieldType dt = initial_dt;
    TimeFieldType t  = current_time();
    assert(Dune::XT::Common::FloatCmp::ge(t_end - t, 0.0));
    size_t time_step_counter = 0;

    const TimeFieldType save_interval = (t_end - t) / num_save_steps;
    TimeFieldType next_save_time      = t + save_interval > t_end ? t_end : t + save_interval;
    size_t save_step_counter          = 1;

    // save/visualize initial solution
    if (save_solution)
      sol.insert(std::make_pair(t, current_solution()));
    if (visualize)
      current_solution().visualize(filename_prefix, Dune::XT::Common::to_string(0));

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      TimeFieldType max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      if (Dune::XT::Common::FloatCmp::gt(t + dt, next_save_time) && num_save_steps != size_t(-1))
        max_dt = std::min(next_save_time - t, max_dt);

      // do a timestep
      dt = step(dt, max_dt);
      t  = current_time();

      // augment time step counter
      ++time_step_counter;

      // check if data should be written in this timestep (and write)
      if (Dune::XT::Common::FloatCmp::ge(t, next_save_time) || num_save_steps == size_t(-1)) {
        if (save_solution)
          sol.insert(sol.end(), std::make_pair(t, current_solution()));
        if (visualize)
          current_solution().visualize(filename_prefix, Dune::XT::Common::to_string(save_step_counter));
        if (output_progress)
          std::cout << "time step " << time_step_counter << " done, time =" << t << ", current dt= " << dt << std::endl;
        next_save_time += save_interval;
        ++save_step_counter;
      }
    } // while (t < t_end)

    return dt;
  } // ... solve(...)

  virtual TimeFieldType solve(const TimeFieldType t_end, const TimeFieldType initial_dt = 1e-4,
                              const size_t num_save_steps = -1, const bool save_solution = true,
                              const bool output_progress = false, const bool visualize = false,
                              const std::string filename_prefix = "solution")
  {
    return solve(
        t_end, initial_dt, num_save_steps, save_solution, output_progress, visualize, filename_prefix, *solution_);
  }

  virtual TimeFieldType solve(const TimeFieldType t_end, const TimeFieldType initial_dt, const size_t num_save_steps,
                              SolutionType& sol)
  {
    return solve(t_end, initial_dt, num_save_steps, true, false, false, "", sol);
  }

  /**
   * \brief Find discrete solution for time point that is closest to t.
   *
   * The timestepper only stores the solution at discrete time points. This function returns the discrete solution for
   * the stored time point that is closest to the query time t.
   *
   * \param t Time
   * \return std::pair with pair.second the discrete solution at time pair.first
   */
  virtual const TimeAndDiscreteFunctionPairType& solution_closest_to_time(const TimeFieldType t) const
  {
    if (solution().empty())
      DUNE_THROW(Dune::InvalidStateException, "Solution is empty!");
    return solution().upper_bound(t)->first - t > t - solution().lower_bound(t)->first ? *solution().lower_bound(t)
                                                                                       : *solution().upper_bound(t);
  }

  virtual const DiscreteFunctionType& solution_at_time(const TimeFieldType t) const
  {
    const auto it = solution().find(t);
    if (it == solution().end())
      DUNE_THROW(Dune::InvalidStateException,
                 "There is no solution for time " + Dune::XT::Common::to_string(t) + " stored in this object!");
    return it->second;
  }

  template <size_t factor = 0>
  void visualize_factor_of_solution(const std::string prefix = "") const
  {
    size_t counter = 0;
    for (const auto& pair : solution()) {
      pair.second.template visualize_factor<factor>(
          prefix + "factor_" + Dune::XT::Common::to_string(factor), Dune::XT::Common::to_string(counter), true);
      ++counter;
    }
  }

  virtual void visualize_solution(const std::string prefix = "") const
  {
    size_t counter = 0;
    for (const auto& pair : solution()) {
      pair.second.visualize(prefix, Dune::XT::Common::to_string(counter));
      ++counter;
    }
  }

private:
  TimeFieldType t_;
  DiscreteFunctionType* u_n_;
  SolutionType* solution_;
}; // class TimeStepperInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_INTERFACE_HH
