// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_TIMESTEPPER_INTERFACE_HH
#define DUNE_GDT_TIMESTEPPER_INTERFACE_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/la/container.hh>


namespace Dune {
namespace GDT {
namespace TimeStepper {

namespace internal {


struct FloatCmpLt
{
  bool operator()(const double& a, const double& b) const
  {
    return DSC::FloatCmp::lt(a, b);
  }
};


}


template< class DiscreteFunctionImp, class TimeFieldImp >
class TimeStepperInterface
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef TimeFieldImp        TimeFieldType;

  typedef typename DiscreteFunctionType::DomainFieldType  DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType   RangeFieldType;
  typedef typename std::map< TimeFieldType, DiscreteFunctionType, typename internal::FloatCmpLt > SolutionType;
  typedef typename SolutionType::value_type TimeAndDiscreteFunctionPairType;

  virtual ~TimeStepperInterface() = default;

  /** Perform one time step with time step length dt and return the estimated optimal time step length for the next
   *  time step. For adaptive methods, dt is just the initial time step length, the actual time step taken may be
   *  shorter. To get the actual time step length, in that case, you need to compare current_time() before and after
   *  calling step().
   * */
  virtual TimeFieldType step(const TimeFieldType dt) = 0;

  virtual TimeFieldType current_time() const = 0;

  virtual const DiscreteFunctionType& solution_at_current_time() const = 0;

  virtual void set_solution_at_current_time(const DiscreteFunctionType& discrete_func) = 0;

  const SolutionType& solution() const
  {
    return solution_;
  }

  /**
   * @brief solve Applies time stepping scheme up to time t_end
   * @param initial_dt Initial time step length. Non-adaptive schemes will use this dt for all time steps. Adaptive
   * schemes will use initial_dt as a first guess for the time step length.
   * @param num_save_steps Number of time points that will be stored in solution if save_solution is true and/or
   * visualized if visualize_solution is true. Default is -1, which will store all time steps.
   * @param save_solution If true, the discrete solution at num_save_steps equidistant time points will be stored in
   * solution.
   * @param output_progress If true, the current time and current time step length will be written to std::cout at
   * num_save_steps equidistant time points
   * @param visualize_solution If true, solution will be written to .vtp/.vtu files at num_save_steps equidistant time
   * points
   * @param filename_prefix Prefix of .vtp/.vtu files for visualization
   * @param solution Optional object that solution will be appended to instead of storing it in the time stepper.
   */
  virtual void solve(const TimeFieldType t_end,
                     const TimeFieldType initial_dt,
                     const size_t num_save_steps,
                     const bool save_solution,
                     const bool output_progress,
                     const bool visualize_solution,
                     const std::string filename_prefix,
                     SolutionType& solution)
  {
    TimeFieldType dt = initial_dt;
    TimeFieldType t = current_time();
    assert(DSC::FloatCmp::ge(t_end - t, 0.0));
    size_t time_step_counter = 0;

    const TimeFieldType save_interval = (t_end - t)/num_save_steps;
    TimeFieldType next_save_time = t + save_interval > t_end ? t_end : t + save_interval;
    size_t save_step_counter = 1;

    while (t < t_end - 1e-10)
    {
      // match saving times and t_end exactly
      if (t + dt > next_save_time && num_save_steps != size_t(-1))
        dt = next_save_time - t;
      if (t + dt > t_end)
        dt = t_end - t;

      // do a timestep
      dt = step(dt);
      t = current_time();

      // check if data should be written in this timestep (and write)
      if (DSC::FloatCmp::ge(t, next_save_time) || num_save_steps == size_t(-1)) {
        if (save_solution)
          solution.insert(std::make_pair(t, solution_at_current_time()));
        if (visualize_solution)
          solution_at_current_time().visualize(filename_prefix, DSC::to_string(save_step_counter));
        if (output_progress)
          std::cout << "time step " << time_step_counter << " done, time =" << t << ", current dt= " << dt << std::endl;
        next_save_time += save_interval;
        ++save_step_counter;
      }

      // augment time step counter
      ++time_step_counter;
    } // while (t < t_end)
  } // ... solve(...)

  virtual void solve(const TimeFieldType t_end,
                     const TimeFieldType initial_dt = 1e-4,
                     const size_t num_save_steps = -1,
                     const bool save_solution = true,
                     const bool output_progress = false,
                     const bool visualize_solution = false,
                     const std::string filename_prefix = "solution")
  {
    solve(t_end, initial_dt, num_save_steps, save_solution, output_progress, visualize_solution, filename_prefix, solution_);
  }

  virtual void solve(const TimeFieldType t_end,
                     const TimeFieldType initial_dt,
                     const size_t num_save_steps,
                     SolutionType& solution)
  {
    solve(t_end, initial_dt, num_save_steps, true, false, false, "", solution);
  }

  virtual const TimeAndDiscreteFunctionPairType& solution_closest_to_time(const TimeFieldType t) const
  {
    if (solution().empty())
      DUNE_THROW(Dune::InvalidStateException, "Solution is empty!");
    return solution().upper_bound(t)->first - t > t - solution().lower_bound(t)->first ? *solution().lower_bound(t) : *solution().upper_bound(t);
  }

  virtual const DiscreteFunctionType& solution_at_time(const TimeFieldType t) const
  {
    const auto it = solution().find(t);
    if (it == solution().end())
        DUNE_THROW(Dune::InvalidStateException, "There is no solution for time " + DSC::to_string(t) + " stored in this object!");
    return it->second;
  }

  template< size_t factor = 0 >
  void visualize_factor_of_solution(const std::string prefix = "") const
  {
    size_t counter = 0;
    for (const auto& pair : solution()) {
      pair.second.template visualize_factor< factor >(prefix + "factor_" + DSC::to_string(factor), DSC::to_string(counter), true);
      ++counter;
    }
  }

  virtual void visualize_solution(const std::string prefix = "") const
  {
    size_t counter = 0;
    for (const auto& pair : solution()) {
      pair.second.visualize(prefix, DSC::to_string(counter));
      ++counter;
    }
  }

private:
  SolutionType solution_;
};

} // namespace TimeStepper
} // namespace Stuff
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_INTERFACE_HH
