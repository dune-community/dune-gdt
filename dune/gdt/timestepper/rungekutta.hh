// This file is part of the dune-stuff project:
//   https://github.com/wwu-numerik/dune-stuff/
// Copyright holders: Rene Milk, Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_STUFF_PLAYGROUND_TIMESTEPPER_HH
#define DUNE_STUFF_PLAYGROUND_TIMESTEPPER_HH

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/string.hh>
#include <dune/stuff/la/container.hh>


namespace Dune {
namespace Stuff {

template< class OperatorImp, class DiscreteFunctionImp, class SourceFunctionImp, size_t stages >
class RungeKuttaTimeStepper
{
  typedef OperatorImp OperatorType;
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef SourceFunctionImp SourceFunctionType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

public:

  RungeKuttaTimeStepper(OperatorType& space_operator,
                        const Dune::FieldMatrix< RangeFieldType, stages + 1, stages + 1 >& butcher_array,
                        const DiscreteFunctionType& initial_values,
                        const SourceFunctionType& source_function)
    : space_operator_(space_operator)
    , butcher_array_(butcher_array)
    , initial_values_(initial_values)
    , u_(initial_values_)
    , source_function_(source_function)
  {
    // copy u_ stages times
    for (size_t ii = 0; ii < stages; ++ii) {
      u_intermediate_stages_.emplace_back(u_);
    }
  }

  void step(const double t_end, const double dt, const double save_interval, const double t_0 = 0)
  {
    double t = t_0;
    int time_step_counter = 0;
    double save_step = save_interval;
    int save_step_counter = 1;

    std::cout << "Visualizing initial values..." << std::endl;
    u_.visualize("concentration_0", false);

    DiscreteFunctionType u_tmp(u_);
    DiscreteFunctionType discrete_source(u_);

    while (t < t_end)
    {
      // evaluate conservation law u_t + L(q) = 0
      for (size_t ii = 0; ii < stages; ++ii) {
        u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
        u_tmp.vector() = u_.vector();
        for (size_t jj = 0; jj < stages; ++jj) {
          u_tmp.vector() += u_intermediate_stages_[jj].vector()*dt*butcher_array_[ii][jj + 1];
        }
        space_operator_.apply(u_tmp , u_intermediate_stages_[ii]);
      };

      for (size_t ii = 0; ii < stages; ++ii) {
        u_.vector() += u_intermediate_stages_[ii].vector()*(-1.0*dt)*butcher_array_[stages][ii+1];
      }

      // evaluate source terms u_t = q(u)
      for (size_t ii = 0; ii < discrete_source.vector().size(); ++ii) {
        discrete_source.vector()[ii] = source_function_.evaluate(u_.vector()[ii]);
      }
      u_.vector() += discrete_source.vector()*dt;

      // augment time step counter
      ++time_step_counter;

      // augment time
      t += dt;

      // check if data should be written
      if (t >= save_step)
      {
        // write data
        u_.visualize("concentration_" + DSC::toString(save_step_counter), false);

        // increase counter and saveStep for next interval
        save_step += save_interval;
        ++save_step_counter;
      }

      // print info about time, timestep size and counter
      std::cout << " k=" << time_step_counter << " t=" << t << " dt=" << dt << std::endl;
    } // while (t < t_end)
  }

private:
  OperatorType& space_operator_;
  const Dune::FieldMatrix< RangeFieldType, stages + 1, stages + 1 >& butcher_array_;
  const DiscreteFunctionType& initial_values_;
  DiscreteFunctionType u_;
  std::vector< DiscreteFunctionType > u_intermediate_stages_;
  const SourceFunctionType& source_function_;
};

} // namespace Stuff
} // namespace Dune

#endif // DUNE_STUFF_PLAYGROUND_TIMESTEPPER_HH
