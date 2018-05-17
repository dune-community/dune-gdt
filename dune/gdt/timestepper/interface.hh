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

#ifndef DUNE_GDT_TIMESTEPPER_INTERFACE_HH
#define DUNE_GDT_TIMESTEPPER_INTERFACE_HH

#include <utility>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/tuple.hh>

#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/datahandle.hh>
#include <dune/gdt/operators/interfaces.hh>

#include "enums.hh"

namespace Dune {
namespace GDT {
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


template <class DiscreteFunctionImp>
class TimeStepperInterface : Dune::XT::Common::StorageProvider<DiscreteFunctionImp>,
                             Dune::XT::Common::StorageProvider<std::map<typename DiscreteFunctionImp::RangeFieldType,
                                                                        DiscreteFunctionImp,
                                                                        typename internal::FloatCmpLt>>
{
public:
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::RangeType RangeType;
  typedef typename std::map<RangeFieldType, DiscreteFunctionType, typename internal::FloatCmpLt> SolutionType;
  typedef typename SolutionType::value_type TimeAndDiscreteFunctionPairType;
  typedef DiscreteFunctionDataHandle<DiscreteFunctionType> DataHandleType;
  typedef std::function<void(const DiscreteFunctionType&, const std::string&, const size_t)> VisualizerType;
  typedef std::function<std::string(const RangeType&)> StringifierType;
  using ThisType = TimeStepperInterface;

private:
  typedef typename Dune::XT::Common::StorageProvider<DiscreteFunctionImp> CurrentSolutionStorageProviderType;
  typedef typename Dune::XT::Common::
      StorageProvider<std::map<RangeFieldType, DiscreteFunctionImp, typename internal::FloatCmpLt>>
          SolutionStorageProviderType;

protected:
  TimeStepperInterface(const RangeFieldType t_0, DiscreteFunctionType& initial_values)
    : CurrentSolutionStorageProviderType(initial_values)
    , SolutionStorageProviderType(new SolutionType())
    , t_(t_0)
    , u_n_(&CurrentSolutionStorageProviderType::access())
    , solution_(&SolutionStorageProviderType::access())
  {
  }

public:
  TimeStepperInterface(const ThisType& other) = delete;
  TimeStepperInterface(ThisType&& source) = delete;
  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

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
  virtual RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) = 0;

  const RangeFieldType& current_time() const
  {
    return t_;
  }

  RangeFieldType& current_time()
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
  virtual RangeFieldType solve(const RangeFieldType t_end,
                               const RangeFieldType initial_dt,
                               const size_t num_save_steps,
                               const bool save_solution,
                               const bool output_progress,
                               const bool visualize,
                               const bool write_to_file,
                               const std::string filename_prefix,
                               SolutionType& sol,
                               const VisualizerType& visualizer,
                               const StringifierType& stringifier)
  {
    RangeFieldType dt = initial_dt;
    RangeFieldType t = current_time();
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    size_t time_step_counter = 0;

    const RangeFieldType save_interval = (t_end - t) / num_save_steps;
    RangeFieldType next_save_time = t + save_interval > t_end ? t_end : t + save_interval;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    if (save_solution)
      sol.insert(std::make_pair(t, current_solution()));
    if (visualize)
      visualizer(current_solution(), filename_prefix, 0);
    if (write_to_file)
      write_to_textfile(current_solution(), filename_prefix, 0, stringifier);

    while (Dune::XT::Common::FloatCmp::lt(t, t_end)) {
      RangeFieldType max_dt = dt;
      // match saving times and t_end exactly
      if (Dune::XT::Common::FloatCmp::gt(t + dt, t_end))
        max_dt = t_end - t;
      if (Dune::XT::Common::FloatCmp::gt(t + dt, next_save_time) && num_save_steps != size_t(-1))
        max_dt = std::min(next_save_time - t, max_dt);

      // do a timestep
      dt = step(dt, max_dt);
      t = current_time();

      // augment time step counter
      ++time_step_counter;

      // check if data should be written in this timestep (and write)
      if (Dune::XT::Common::FloatCmp::ge(t, next_save_time) || num_save_steps == size_t(-1)) {
        if (save_solution)
          sol.insert(sol.end(), std::make_pair(t, current_solution()));
        if (visualize)
          visualizer(current_solution(), filename_prefix, save_step_counter);
        if (write_to_file)
          write_to_textfile(current_solution(), filename_prefix, save_step_counter, stringifier);
        if (output_progress)
          std::cout << "time step " << time_step_counter << " done, time =" << t << ", current dt= " << dt << std::endl;
        next_save_time += save_interval;
        ++save_step_counter;
      }
    } // while (t < t_end)

    return dt;
  } // ... solve(...)

  virtual RangeFieldType solve(const RangeFieldType t_end,
                               const RangeFieldType initial_dt = 1e-4,
                               const size_t num_save_steps = -1,
                               const bool save_solution = true,
                               const bool output_progress = false,
                               const bool visualize = false,
                               const std::string filename_prefix = "solution",
                               const VisualizerType& visualizer = all_components_visualizer(),
                               const StringifierType& stringifier = vector_stringifier())
  {
    return solve(t_end,
                 initial_dt,
                 num_save_steps,
                 save_solution,
                 output_progress,
                 visualize,
                 false,
                 filename_prefix,
                 *solution_,
                 visualizer,
                 stringifier);
  }

  virtual RangeFieldType solve(const RangeFieldType t_end,
                               const RangeFieldType initial_dt,
                               const size_t num_save_steps,
                               const bool save_solution,
                               const bool output_progress,
                               const bool visualize,
                               const bool write_to_file,
                               const std::string filename_prefix = "solution",
                               const VisualizerType& visualizer = all_components_visualizer(),
                               const StringifierType& stringifier = vector_stringifier())
  {
    return solve(t_end,
                 initial_dt,
                 num_save_steps,
                 save_solution,
                 output_progress,
                 visualize,
                 write_to_file,
                 filename_prefix,
                 *solution_,
                 visualizer,
                 stringifier);
  }


  virtual RangeFieldType
  solve(const RangeFieldType t_end, const RangeFieldType initial_dt, const size_t num_save_steps, SolutionType& sol)
  {
    return solve(t_end,
                 initial_dt,
                 num_save_steps,
                 true,
                 false,
                 false,
                 false,
                 "",
                 sol,
                 all_components_visualizer(),
                 vector_stringifier());
  }

  /**
   * \brief Find discrete solution for time point that is closest to t.
   *
   * The timestepper only stores the solution at discrete time points. This function returns the discrete solution
   * for
   * the stored time point that is closest to the query time t.
   *
   * \param t Time
   * \return std::pair with pair.second the discrete solution at time pair.first
   */
  virtual const TimeAndDiscreteFunctionPairType& solution_closest_to_time(const RangeFieldType t) const
  {
    if (solution().empty())
      DUNE_THROW(Dune::InvalidStateException, "Solution is empty!");
    return solution().upper_bound(t)->first - t > t - solution().lower_bound(t)->first ? *solution().lower_bound(t)
                                                                                       : *solution().upper_bound(t);
  }

  virtual const DiscreteFunctionType& solution_at_time(const RangeFieldType t) const
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

  static VisualizerType all_components_visualizer()
  {
    return [](const DiscreteFunctionType& u_n, const std::string& filename_prefix, const size_t step) {
      u_n.visualize(filename_prefix, Dune::XT::Common::to_string(step));
    };
  }

  static StringifierType vector_stringifier()
  {
    return [](const RangeType& val) { return XT::Common::to_string(val, 15); };
  } // ... vector_stringifier()

  void write_to_textfile(const DiscreteFunctionType& u_n,
                         const std::string& filename_prefix,
                         const size_t step,
                         const StringifierType& stringifier)
  {
    const auto& grid_layer = u_n.space().grid_layer();
    std::ofstream valuesfile(filename_prefix + "_values_rank_" + XT::Common::to_string(grid_layer.comm().rank()) + "_"
                             + XT::Common::to_string(step)
                             + ".txt");
    for (const auto& entity : elements(grid_layer, Dune::Partitions::interiorBorder)) {
      const auto local_func = u_n.local_function(entity);
      const auto entity_center = entity.geometry().center();
      const auto val = local_func->evaluate(entity.geometry().local(entity_center));
      for (size_t ii = 0; ii < entity_center.size(); ++ii)
        valuesfile << XT::Common::to_string(entity_center[ii], 15) << " ";
      valuesfile << stringifier(val) << std::endl;
    }
    valuesfile.close();
    grid_layer.comm().barrier();
    if (grid_layer.comm().rank() == 0) {
      std::remove((filename_prefix + "_values_" + XT::Common::to_string(step) + ".txt").c_str());
      std::ofstream merged_valuesfile(filename_prefix + "_values_" + XT::Common::to_string(step) + ".txt",
                                      std::ios_base::binary | std::ios_base::app);
      for (int ii = 0; ii < grid_layer.comm().size(); ++ii) {
        std::ifstream valuefile_rank(
            filename_prefix + "_values_rank_" + XT::Common::to_string(ii) + "_" + XT::Common::to_string(step) + ".txt",
            std::ios_base::binary);
        merged_valuesfile << valuefile_rank.rdbuf();
        valuefile_rank.close();
      } // ii
      merged_valuesfile.close();
    } // if (rank == 0)
  } // void write_to_textfile()

private:
  RangeFieldType t_;
  DiscreteFunctionType* u_n_;
  SolutionType* solution_;
}; // class TimeStepperInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_INTERFACE_HH
