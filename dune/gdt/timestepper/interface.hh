// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TIMESTEPPER_INTERFACE_HH
#define DUNE_GDT_TIMESTEPPER_INTERFACE_HH

#include <utility>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/discretefunction/datahandle.hh>
#include <dune/gdt/operators/interfaces.hh>

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
  explicit_rungekutta_other,
  implicit_euler,
  implicit_midpoint,
  trapezoidal_rule,
  diagonally_implicit_other
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

template <size_t ii, class DiscreteFunctionType>
auto function_factor(const DiscreteFunctionType& discrete_function) -> typename Dune::GDT::DiscreteFunction<
    typename std::tuple_element<ii, typename DiscreteFunctionType::SpaceType::SpaceTupleType>::type,
    typename DiscreteFunctionType::VectorType>
{
  typedef typename Dune::GDT::DiscreteFunction<
      typename std::tuple_element<ii, typename DiscreteFunctionType::SpaceType::SpaceTupleType>::type,
      typename DiscreteFunctionType::VectorType>
      FactorDiscreteFunctionType;
  static_assert(ii < DiscreteFunctionType::SpaceType::num_factors, "This factor does not exist.");
  const auto& space = discrete_function.space();
  const auto& factor_space = space.template factor<ii>();
  typename DiscreteFunctionType::VectorType factor_vector(factor_space.mapper().size());
  const auto it_end = space.grid_view().template end<0>();
  for (auto it = space.grid_view().template begin<0>(); it != it_end; ++it) {
    const auto& entity = *it;
    for (size_t jj = 0; jj < factor_space.mapper().numDofs(entity); ++jj)
      factor_vector[factor_space.mapper().mapToGlobal(entity, jj)] =
          discrete_function.vector()[space.mapper().mapToGlobal(ii, entity, jj)];
  }
  FactorDiscreteFunctionType factor_discrete_function(factor_space);
  factor_discrete_function.vector() = factor_vector;
  typedef Dune::GDT::DiscreteFunctionDataHandle<FactorDiscreteFunctionType> DataHandleType;
  DataHandleType handle(factor_discrete_function);
  factor_space.grid_view().template communicate<DataHandleType>(
      handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
  return factor_discrete_function;
}

template <size_t current_factor_index, size_t last_factor_index>
struct static_for_loop
{
  template <class DiscreteFunctionType, class FactorDiscreteFunctionType>
  static void sum_vector(const DiscreteFunctionType& discrete_function,
                         FactorDiscreteFunctionType& first_discrete_function)
  {
    first_discrete_function.vector() +=
        function_factor<current_factor_index, DiscreteFunctionType>(discrete_function).vector();
    static_for_loop<current_factor_index + 1, last_factor_index>::sum_vector(discrete_function,
                                                                             first_discrete_function);
  }

  template <class DiscreteFunctionType, class FactorDiscreteFunctionType>
  static void sum_even_vector(const DiscreteFunctionType& discrete_function,
                              FactorDiscreteFunctionType& first_discrete_function)
  {
    if (!(current_factor_index % 2))
      first_discrete_function.vector() +=
          function_factor<current_factor_index, DiscreteFunctionType>(discrete_function).vector();
    static_for_loop<current_factor_index + 1, last_factor_index>::sum_even_vector(discrete_function,
                                                                                  first_discrete_function);
  }
};

// specialization of static for loop to end the loop
template <size_t last_factor_index>
struct static_for_loop<last_factor_index, last_factor_index>
{
  template <class DiscreteFunctionType, class FactorDiscreteFunctionType>
  static void sum_vector(const DiscreteFunctionType& discrete_function,
                         FactorDiscreteFunctionType& first_discrete_function)
  {
    first_discrete_function.vector() +=
        function_factor<last_factor_index, DiscreteFunctionType>(discrete_function).vector();
  }

  template <class DiscreteFunctionType, class FactorDiscreteFunctionType>
  static void sum_even_vector(const DiscreteFunctionType& discrete_function,
                              FactorDiscreteFunctionType& first_discrete_function)
  {
    if (!(last_factor_index % 2))
      first_discrete_function.vector() +=
          function_factor<last_factor_index, DiscreteFunctionType>(discrete_function).vector();
  }
};


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
  typedef DiscreteFunctionDataHandle<DiscreteFunctionType> DataHandleType;

private:
  typedef typename Dune::XT::Common::StorageProvider<DiscreteFunctionImp> CurrentSolutionStorageProviderType;
  typedef typename Dune::XT::Common::
      StorageProvider<std::map<TimeFieldImp, DiscreteFunctionImp, typename internal::FloatCmpLt>>
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
  virtual TimeFieldType solve(const TimeFieldType t_end,
                              const TimeFieldType initial_dt,
                              const size_t num_save_steps,
                              const bool save_solution,
                              const bool output_progress,
                              const bool visualize,
                              const std::string filename_prefix,
                              SolutionType& sol,
                              const int visualize_tag = 0)
  {
    TimeFieldType dt = initial_dt;
    TimeFieldType t = current_time();
    assert(Dune::XT::Common::FloatCmp::ge(t_end - t, 0.0));
    size_t time_step_counter = 0;

    const TimeFieldType save_interval = (t_end - t) / num_save_steps;
    TimeFieldType next_save_time = t + save_interval > t_end ? t_end : t + save_interval;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    if (save_solution)
      sol.insert(std::make_pair(t, current_solution()));
    if (visualize) {
      const auto dimRange = DiscreteFunctionType::dimRange;
      const auto& u_n = current_solution();
      if (visualize_tag == 0) {
        u_n.visualize(filename_prefix, Dune::XT::Common::to_string(0));
      } else if (visualize_tag == 1) {
        auto sum_function = function_factor<0, DiscreteFunctionType>(u_n);
        static_for_loop<1, dimRange - 1>::sum_vector(u_n, sum_function);
        sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(0));
      } else if (visualize_tag == 2) {
        auto sum_function = function_factor<0, DiscreteFunctionType>(u_n);
        static_for_loop<1, dimRange - 1>::sum_even_vector(u_n, sum_function);
        sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(0));
      }
    }


    while (Dune::XT::Common::FloatCmp::lt(t, t_end, 1e-10)) {
      TimeFieldType max_dt = dt;
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
        if (visualize) {
          const auto dimRange = DiscreteFunctionType::dimRange;
          const auto& u_n = current_solution();
          if (visualize_tag == 0) {
            u_n.visualize(filename_prefix, Dune::XT::Common::to_string(save_step_counter));
          } else if (visualize_tag == 1) {
            auto sum_function = function_factor<0, DiscreteFunctionType>(u_n);
            static_for_loop<1, dimRange - 1>::sum_vector(u_n, sum_function);
            sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(save_step_counter));
          } else if (visualize_tag == 2) {
            auto sum_function = function_factor<0, DiscreteFunctionType>(u_n);
            static_for_loop<1, dimRange - 1>::sum_even_vector(u_n, sum_function);
            sum_function.visualize(filename_prefix + "_" + Dune::XT::Common::to_string(save_step_counter));
          }
        }
        if (output_progress)
          std::cout << "time step " << time_step_counter << " done, time =" << t << ", current dt= " << dt << std::endl;
        next_save_time += save_interval;
        ++save_step_counter;
      }
    } // while (t < t_end)

    return dt;
  } // ... solve(...)

  virtual TimeFieldType solve(const TimeFieldType t_end,
                              const TimeFieldType initial_dt = 1e-4,
                              const size_t num_save_steps = -1,
                              const bool save_solution = true,
                              const bool output_progress = false,
                              const bool visualize = false,
                              const std::string filename_prefix = "solution",
                              const int visualize_prefix = 0)
  {
    return solve(t_end,
                 initial_dt,
                 num_save_steps,
                 save_solution,
                 output_progress,
                 visualize,
                 filename_prefix,
                 *solution_,
                 visualize_prefix);
  }

  virtual TimeFieldType
  solve(const TimeFieldType t_end, const TimeFieldType initial_dt, const size_t num_save_steps, SolutionType& sol)
  {
    return solve(t_end, initial_dt, num_save_steps, true, false, false, "", sol, 0);
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
