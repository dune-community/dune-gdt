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
 * and q an operator representing a source or sink.
 * A fractional step approach is used to evolve the equation, where the same Runge Kutta method is used in both steps.
 * The specific Runge Kutta method can be chosen in the constructor by supplying a DynamicMatrix< RangeFieldType >
 * A and vectors (DynamicVector< RangeFieldType >) b and c. Here, A, b and c form the butcher tableau (see
 * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods, A is composed of the coefficients a_{ij}, b of b_j
 * and c of c_j). The default is a forward euler method. By now, c will be ignored as operators that are explicitly
 * time-dependent are not supported yet.
 *
 * \tparam FluxOperatorImp Type of flux operator L
 * \tparam DiscreteFunctionImp Type of initial values
 * \tparam SourceFunctionImp Type of source operator q
 */

template< class FluxOperatorImp, class SourceOperatorImp, class DiscreteFunctionImp, class TimeFieldImp >
class RungeKutta
{
public:
  typedef FluxOperatorImp     FluxOperatorType;
  typedef SourceOperatorImp   SourceOperatorType;
  typedef DiscreteFunctionImp DiscreteFunctionType;
  typedef TimeFieldImp        TimeFieldType;

  typedef typename DiscreteFunctionType::DomainFieldType  DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType   RangeFieldType;
  typedef typename Dune::DynamicMatrix< RangeFieldType >  MatrixType;
  typedef typename Dune::DynamicVector< RangeFieldType >  VectorType;
  typedef typename Dune::DynamicVector< TimeFieldType >   TimeVectorType;
  typedef typename std::vector< std::pair< TimeFieldType, DiscreteFunctionType > > SolutionType;

  /**
   * \brief Constructor for RungeKutta time stepper
   *
   * \param flux operator L
   * \param source operator q
   * \param initial_values Discrete function containing initial values for u
   * \param start_time Starting time (s.t. u(start_time) = initial_values)
   * \param A A (see above)
   * \param b b (see above)
   * \param c c (completely ignored, see above)
   */
  RungeKutta(const FluxOperatorType& flux_operator,
             const SourceOperatorType& source_operator,
             const DiscreteFunctionType& initial_values,
             const DomainFieldType dx,
             const MatrixType A = DSC::fromString< MatrixType >("[0]"),
             const VectorType b = DSC::fromString< VectorType >("[1]"),
             const TimeVectorType c = DSC::fromString< TimeVectorType >("0"))
    : flux_operator_(flux_operator)
    , source_operator_(source_operator)
    , initial_values_(initial_values)
    , u_n_(initial_values_)
    , u_tmp_(u_n_)
    , t_(0.0)
    , dx_(dx)
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
                 "A has to be a lower triangular matrix with 0 on the diagonal (implicit methods are not implemented)");
        }
    }
#endif //NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_ ; ++ii) {
      u_intermediate_stages_.emplace_back(u_n_);
    }
  } // constructor

  TimeFieldType step(const TimeFieldType dt)
  {
//    DiscreteFunctionType u_n_copy(u_n_);
    apply_RK_scheme(flux_operator_, dt, -1.0);       // evaluate conservation law d_t u + L(u) = 0
    apply_RK_scheme(source_operator_,dt, 1.0);      // evaluate source terms d_t u = q(u)

      // unsplit method (comment second apply_RK_scheme above and uncomment definition of u_n_copy)
//      u_intermediate_stages_[0].vector() *= RangeFieldType(0);
//      source_operator_.apply(u_n_copy, u_intermediate_stages_[0], t_);
//      u_n_.vector() += u_intermediate_stages_[0].vector()*dt;

      t_ += dt;                                        // augment time

      //calculate new dt <= dx/(2*max_j abs(u_j)) (for TVD MUSCL, see FiniteVolumenLiteratur/TVD-RungeKutta-Schemes)
//      RangeFieldType max_u_j_abs = 0;
//      for (auto& u_j : u_n_.vector()) {
//        const RangeFieldType u_j_abs = std::abs(u_j);
//        if (u_j_abs > max_u_j_abs)
//          max_u_j_abs = u_j_abs;
//      }
      TimeFieldType dt_new = dt; //0.99*dx_/(8.0*max_u_j_abs);

      // return
      return dt_new;
  } // ... step(...)

  template< class OperatorImp >
  void apply_RK_scheme(const OperatorImp& op, const TimeFieldType dt, const RangeFieldType factor)
  {
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_intermediate_stages_[ii].vector() *= RangeFieldType(0);
      u_tmp_.vector() = u_n_.vector();
      for (size_t jj = 0; jj < ii; ++jj)
        u_tmp_.vector() += u_intermediate_stages_[jj].vector()*(dt*factor*(A_[ii][jj]));
      op.apply(u_tmp_, u_intermediate_stages_[ii], t_ + dt*c_[ii]);
    }

    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_n_.vector() += u_intermediate_stages_[ii].vector()*(factor*dt*b_[ii]);
    }
  } // void apply_RK_scheme(...)

  void reset()
  {
    t_ = 0.0;
    u_n_.vector() = initial_values_.vector();
  }

  void solve(const TimeFieldType t_end,
             const TimeFieldType first_dt,
             const TimeFieldType save_step_length,
             const bool save_solution,
             const bool write_solution,
             const std::string filename_prefix,
             std::vector< std::pair< double, DiscreteFunctionType > >& solution)
  {
    TimeFieldType dt = first_dt;
    assert(t_end - t_ >= dt);
    size_t time_step_counter = 0;

    const TimeFieldType save_interval = DSC::FloatCmp::eq(save_step_length, 0.0) ? dt : save_step_length;
    const TimeFieldType output_interval = 0.001;
    TimeFieldType next_save_time = t_ + save_interval > t_end ? t_end : t_ + save_interval;
    TimeFieldType next_output_time = t_ + output_interval;
    size_t save_step_counter = 1;

    // clear solution
    if (save_solution) {
      solution.clear();
      solution.emplace_back(std::make_pair(t_, u_n_));
    }
    if (write_solution)
      u_n_.template visualize_factor< 0 >(filename_prefix + "factor_0_0");

    while (t_ + dt < t_end)
    {
      // do a timestep
      dt = step(dt);

      // check if data should be written in this timestep (and write)
      if (DSC::FloatCmp::ge(t_, next_save_time - 1e-10)) {
        if (save_solution)
          solution.emplace_back(std::make_pair(t_, u_n_));
        if (write_solution)
          u_n_.template visualize_factor< 0 >(filename_prefix + "factor_0_" + DSC::toString(save_step_counter));
        next_save_time += save_interval;
        ++save_step_counter;
      }

      // augment time step counter
      ++time_step_counter;

      // print info about time, timestep size and counter
//      if (DSC::FloatCmp::ge(t_, next_output_time)) {
//        std::cout << " k=" << time_step_counter << " t=" << t_ << " dt=" << dt << std::endl;
//        next_output_time += output_interval;
//      }
    } // while (t_ < t_end)

    // do last step s.t. it matches t_end exactly
    if (!DSC::FloatCmp::ge(t_, t_end - 1e-10)) {
      step(t_end - t_);
      solution.emplace_back(std::make_pair(t_, u_n_));
    }
  } // ... solve(...)

  void solve(const TimeFieldType t_end,
             const TimeFieldType first_dt,
             const TimeFieldType save_step_length = 0.0,
             const bool save_solution = false,
             const bool write_solution = true,
             const std::string filename_prefix = "")
  {
    solve(t_end, first_dt, save_step_length, save_solution, write_solution, filename_prefix, solution_);
  }

  void solve(const TimeFieldType t_end,
             const TimeFieldType first_dt,
             const TimeFieldType save_step_length,
             std::vector< std::pair< double, DiscreteFunctionType > >& solution)
  {
    solve(t_end, first_dt, save_step_length, true, false, "", solution);
  }

  TimeFieldType current_time() const
  {
    return t_;
  }

  const DiscreteFunctionType& current_solution() const
  {
    return u_n_;
  }

  template< size_t factor_to_be_visualized = 0 >
  void visualize_solution(const std::string prefix = "") const
  {
    for (size_t ii = 0; ii < solution_.size(); ++ii) {
      auto& pair = solution_[ii];
      pair.second.template visualize_factor< factor_to_be_visualized >(prefix + "factor_"
                                                                       + DSC::toString(factor_to_be_visualized)
                                                                       + "_" + DSC::toString(ii), true);
    }
  }

  const std::pair< bool, TimeFieldType > find_suitable_dt(const TimeFieldType initial_dt,
                                                          const TimeFieldType dt_refinement_factor = 2,
                                                          const RangeFieldType treshold = 0.9*std::numeric_limits< RangeFieldType >::max(),
                                                          const size_t max_steps_per_dt = 20,
                                                          const size_t max_refinements = 20)
  {
    assert(treshold > 0);
    // save current state
    DiscreteFunctionType initial_u_n = u_n_;
    TimeFieldType initial_t = t_;
    // start with initial dt
    TimeFieldType current_dt = initial_dt;
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

  const SolutionType solution() const {
    return solution_;
  }

private:
  const FluxOperatorType& flux_operator_;
  const SourceOperatorType& source_operator_;
  const DiscreteFunctionType& initial_values_;
  DiscreteFunctionType u_n_;
  DiscreteFunctionType u_tmp_;
  double t_;
  double dx_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector< DiscreteFunctionType > u_intermediate_stages_;
  const size_t num_stages_;
  SolutionType solution_;
};

} // namespace TimeStepper
} // namespace Stuff
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_RUNGEKUTTA_HH
