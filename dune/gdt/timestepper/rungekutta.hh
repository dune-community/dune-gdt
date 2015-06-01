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


namespace Dune {
namespace GDT {
namespace TimeStepper {

/** \brief Time stepper using Runge Kutta methods
 *
 * Timestepper for equations of the form u_t + L(u) = q(u) where u is a discrete function, L a space operator
 * and q a function representing a source or sink.
 * A fractional step approach is used to evolve the equation, where the same Runge Kutta method is used in both steps.
 * The specific Runge Kutta method can be chosen in the constructor by supplying a DynamicMatrix< RangeFieldType >
 * A and vectors (DynamicVector< RangeFieldType >) b and c. Here, A, b and c form the butcher tableau (see
 * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods, A is composed of the coefficients a_{ij}, b of b_j
 * and c of c_j). The default is a forward euler method. By now, c will be ignored as operators that are explicitly
 * time-dependent are not supported yet.
 *
 * \tparam OperatorImp Type of space operator, has to offer a void apply(DiscreteFunctionImp, std::vector) method
 * \tparam DiscreteFunctionImp Type of initial values
 * \tparam SourceFunctionImp Type of source function, has to offer a RangeType evaluate(DomainType) method
 */
template< class OperatorImp, class DiscreteFunctionImp, class SourceFunctionImp >
class RungeKutta
{
  typedef OperatorImp OperatorType;
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef SourceFunctionImp SourceFunctionType;

public:
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::DynamicMatrix< RangeFieldType > MatrixType;
  typedef typename Dune::DynamicVector< RangeFieldType > VectorType;

  /**
   * \brief Constructor for RungeKutta time stepper
   *
   * \param space_operator L
   * \param initial_values Discrete function containing initial values for u
   * \param source_function q
   * \param start_time Starting time (s.t. u(start_time) = initial_values)
   * \param A A (see above)
   * \param b b (see above)
   * \param c c (completely ignored, see above)
   */
  RungeKutta(OperatorType& space_operator,
             const DiscreteFunctionType& initial_values,
             const SourceFunctionType& source_function,
             const MatrixType A = DSC::fromString< MatrixType >("[0]"),
             const VectorType b = DSC::fromString< VectorType >("[1]"),
             const VectorType c = DSC::fromString< VectorType >("[0]"))
    : space_operator_(space_operator)
    , initial_values_(initial_values)
    , u_n_(initial_values_)
    , source_function_(source_function)
    , t_(0.0)
    , A_(A)
    , b_(b)
    , c_(c)
    , num_stages_(A_.rows())
    , solution_(std::vector< std::pair< double, DiscreteFunctionType > >())
  {
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_.size() == A_.rows());
//    assert(c_.size() == A_.rows());
#ifndef NDEBUG
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
        for (size_t jj = ii; jj < A_.cols(); ++jj) {
          assert(DSC::FloatCmp::eq(A_[ii][jj], 0.0) &&
                 "A has to be a lower triangular matrix with 0 on the diagonal (implicit methods are not supported)");
        }
    }
#endif //NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_ ; ++ii) {
      u_intermediate_stages_.emplace_back(u_n_);
    }
  } // constructor

  double step(const double dt)
  {
      DiscreteFunctionType u_tmp(u_n_);

      // evaluate conservation law u_t + L(u) = 0
      for (size_t ii = 0; ii < num_stages_; ++ii) {
        u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
        u_tmp.vector() = u_n_.vector();
        for (size_t jj = 0; jj < num_stages_; ++jj) {
          u_tmp.vector() += u_intermediate_stages_[jj].vector()*dt*A_[ii][jj];
        }
        space_operator_.apply(u_tmp , u_intermediate_stages_[ii], t_);
      };

      for (size_t ii = 0; ii < num_stages_; ++ii) {
        u_n_.vector() += u_intermediate_stages_[ii].vector()*(-1.0*dt)*b_[ii];
      }

      // evaluate source terms u_t = q(u)
      for (size_t ii = 0; ii < num_stages_; ++ii) {
        u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
        u_tmp.vector() = u_n_.vector();
        for (size_t jj = 0; jj < num_stages_; ++jj) {
          u_tmp.vector() += u_intermediate_stages_[jj].vector()*dt*A_[ii][jj];
        }
        const auto it_end = u_n_.space().grid_view().template end< 0 >();
        for (auto it = u_n_.space().grid_view().template begin< 0 >(); it != it_end; ++it) {
          const auto& entity = *it;
          const auto source_value = source_function_.local_global_function(entity)->evaluate(entity.geometry().local(entity.geometry().center()), u_tmp.local_function(entity)->evaluate(entity.geometry().local(entity.geometry().center())));
          for (size_t kk = 0; kk < source_value.size(); ++kk)
            u_intermediate_stages_[ii].local_discrete_function(entity)->vector().set(kk, source_value[kk]);
        }
      };

      for (size_t ii = 0; ii < num_stages_; ++ii) {
        u_n_.vector() += u_intermediate_stages_[ii].vector()*dt*b_[ii];
      }

      // augment time
      t_ += dt;

      // return
      return dt;
  } // ... step(...)

  void solve(const double t_end,
             const double first_dt,
             const double save_step,
             const bool output,
             const bool save_solution,
             std::vector< std::pair< double, DiscreteFunctionType > >& solution)
  {
    double dt = first_dt;
    assert(t_end - t_ >= dt);
    size_t time_step_counter = 0;

    const double save_interval = DSC::FloatCmp::eq(save_step, 0.0) ? dt : save_step;
    double next_save_time = t_ + save_interval;
    size_t save_step_counter = 1;

    const size_t factor_to_be_visualized = 0;

    if (output)
      std::cout << "Visualizing initial values..." << std::endl;
    u_n_.template visualize_factor< factor_to_be_visualized >("factor_" + DSC::toString(factor_to_be_visualized) + "_0", false);

    if (save_solution) {
      // clear solution
      solution.clear();
      solution.emplace_back(std::make_pair(t_, u_n_));
    }

    if (output)
      std::cout << "Starting time loop..." << std::endl;
    while (t_ < t_end)
    {
      // do a timestep
      dt = step(dt);

      // check if data should be written in this timestep (and write)
      if (t_ >= next_save_time) {
        u_n_.template visualize_factor< factor_to_be_visualized >("factor_" + DSC::toString(factor_to_be_visualized) + "_" + DSC::toString(save_step_counter), false);
        if (save_solution) {
          solution.emplace_back(std::make_pair(t_, u_n_));
        }
        next_save_time += save_interval;
        ++save_step_counter;
      }

      // augment time step counter
      ++time_step_counter;

      // print info about time, timestep size and counter
      if (output)
        std::cout << " k=" << time_step_counter << " t=" << t_ << " dt=" << dt << std::endl;
    } // while (t_ < t_end)
  } // ... solve(...)

  void solve(const double t_end,
             const double first_dt,
             const double save_step = 0.0,
             const bool output = false,
             const bool save_solution = false)
  {
    solve(t_end, first_dt, save_step, output, save_solution, solution_);
  }

  double current_time() const
  {
    return t_;
  }

  const DiscreteFunctionType& current_solution() const
  {
    return u_n_;
  }

  const std::pair< bool, double > find_suitable_dt(const double initial_dt,
                                                   const double dt_refinement_factor = 2,
                                                   const double treshold = 0.9*std::numeric_limits< double >::max(),
                                                   const size_t max_steps_per_dt = 20,
                                                   const size_t max_refinements = 20)
  {
    assert(treshold > 0);
    // save current state
    DiscreteFunctionType initial_u_n = u_n_;
    double initial_t = t_;
    // start with initial dt
    double current_dt = initial_dt;
    size_t num_refinements = 0;
    while (num_refinements < max_refinements) {
      std::cout << "Trying time step length dt = " << current_dt << "... ";
      bool unlikely_value_occured = false;
      size_t num_steps = 0;
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
  OperatorType& space_operator_;
  const DiscreteFunctionType& initial_values_;
  DiscreteFunctionType u_n_;
  const SourceFunctionType& source_function_;
  double t_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector< DiscreteFunctionType > u_intermediate_stages_;
  const size_t num_stages_;
  std::vector< std::pair< double, DiscreteFunctionType > > solution_;
};

} // namespace TimeStepper
} // namespace Stuff
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_RUNGEKUTTA_HH
