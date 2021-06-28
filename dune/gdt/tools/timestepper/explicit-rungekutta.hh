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
    return {{0.}};
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return {1.};
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return {0.};
  }
};

// Second order SSP
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_second_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return {{0., 0.}, {1., 0.}};
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return {0.5, 0.5};
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return {0., 1.};
  }
};

// Third order SSP
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_third_order_ssp>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return {{0., 0., 0.}, {1., 0., 0.}, {0.25, 0.25, 0.}};
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return {1. / 6., 1. / 6., 2. / 3.};
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return {0., 1., 0.5};
  }
};

// Classic fourth order RK
template <class RangeFieldType>
struct ButcherArrayProvider<RangeFieldType, TimeStepperMethods::explicit_rungekutta_classic_fourth_order>
{
  static Dune::DynamicMatrix<RangeFieldType> A()
  {
    return {{0., 0., 0., 0.}, {0.5, 0., 0., 0.}, {0., 0.5, 0., 0.}, {0., 0., 1., 0}};
  }

  static Dune::DynamicVector<RangeFieldType> b()
  {
    return {1. / 6., 1. / 3., 1. / 3., 1. / 6.};
  }

  static Dune::DynamicVector<RangeFieldType> c()
  {
    return {0., 0.5, 0.5, 1.};
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
  using BaseType::dimRange;

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
    , op_(op)
    , r_(r)
    , u_i_(BaseType::current_solution().copy_as_discrete_function())
    , A_(A)
    , b_(b)
    , c_(c)
    , num_stages_(A_.rows())
    , relaxationupdate_stages_(num_stages_)
    , last_entropy_(0)
  {
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_.size() == A_.rows());
    assert(c_.size() == A_.rows());
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        DUNE_THROW_IF(XT::Common::FloatCmp::ne(A_[ii][jj], 0.0),
                      XT::Common::Exceptions::wrong_input_given,
                      "A has to be a lower triangular matrix with 0 on the main diagonal");
      }
    }
    // store as many discrete functions as needed for the stages k
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(current_solution().copy_as_discrete_function());
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

  RangeFieldType
  step(const RangeFieldType dt, const RangeFieldType max_dt, const std::string prefix = "") override final
  {
    auto& u_n = current_solution();
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    this->dts_.push_back(actual_dt);

    const auto local_u = u_n.local_discrete_function();
    std::vector<double> val_vector;
    const auto& grid_view = u_n.space().grid_view();
    const auto& quadratures = op_.entropy_solver().entropy_flux().basis_functions().quadratures();
    static const auto merged_quads = XT::Data::merged_quadrature(quadratures);
    static const auto basis_vals = get_basis_vals(merged_quads);
    if (XT::Common::FloatCmp::eq(t, 0.)) {
      last_entropy_ = compute_entropy(local_u, grid_view, merged_quads, basis_vals, val_vector);
      write_entropy(local_u, grid_view, 0., merged_quads, basis_vals, t, actual_dt, val_vector, prefix);
    }

    // calculate stages
    double relaxationupdate = 0.;
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      u_i_->dofs().vector() = u_n.dofs().vector();
      for (size_t jj = 0; jj < ii; ++jj)
        u_i_->dofs().vector() += stages_k_[jj]->dofs().vector() * (actual_dt * r_ * (A_[ii][jj]));
      // TODO: provide actual_dt to op_. This leads to spurious oscillations in the Lax-Friedrichs flux
      // because actual_dt/dx may become very small.
      relaxationupdate_stages_[ii] = 0.;
      op_.apply_and_compute_relax(u_i_->dofs().vector(),
                                  stages_k_[ii]->dofs().vector(),
                                  XT::Common::Parameter({{"t", {t + actual_dt * c_[ii]}}, {"dt", {dt}}}),
                                  relaxationupdate_stages_[ii]);
      relaxationupdate += b_[ii] * relaxationupdate_stages_[ii];
      DataHandleType stages_k_ii_handle(*stages_k_[ii]);
      stages_k_[ii]->space().grid_view().template communicate<DataHandleType>(
          stages_k_ii_handle, Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    }

    // calculate value of u at next time step
    for (size_t ii = 0; ii < num_stages_; ++ii)
      u_n.dofs().vector() += stages_k_[ii]->dofs().vector() * (r_ * actual_dt * b_[ii]);

    // augment time
    t += actual_dt;

    write_entropy(local_u, grid_view, relaxationupdate, merged_quads, basis_vals, t, actual_dt, val_vector, prefix);

    return dt;
  } // ... step(...)

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
  template <class MergedQuadratureType>
  std::vector<XT::Common::FieldVector<double, dimRange>> get_basis_vals(const MergedQuadratureType& merged_quads)
  {
    const auto& basis_functions = op_.entropy_solver().entropy_flux().basis_functions();
    std::vector<XT::Common::FieldVector<double, dimRange>> basis_vals(merged_quads.size());
    size_t kk = 0;
    for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it, ++kk) {
      const auto& quad_point = *it;
      const auto& v = quad_point.position();
      basis_vals[kk] = basis_functions.evaluate(v, it.first_index());
    }
    return basis_vals;
  }

  template <class GV, class LocalDiscreteFunctionType, class MergedQuadratureType, class BasisValsType>
  void write_entropy(LocalDiscreteFunctionType& local_u,
                     const GV& grid_view,
                     const double relaxationupdate,
                     const MergedQuadratureType& merged_quads,
                     const BasisValsType& basis_vals,
                     const double t,
                     const double dt,
                     std::vector<double>& val_vector,
                     std::string prefix = "entropy")
  {
    const std::string entropy_filename = prefix + "_entropy.txt";
    std::ofstream entropy_file(entropy_filename, std::ios_base::app);
    const double entropy = compute_entropy(local_u, grid_view, merged_quads, basis_vals, val_vector);
    const double diff =
        XT::Common::FloatCmp::eq(t, 0.) ? 0. : compensated_sum({entropy, -1. * last_entropy_, -dt * relaxationupdate});
    entropy_file << XT::Common::to_string(t) << " " << XT::Common::to_string(entropy, 15) << " "
                 << XT::Common::to_string(diff, 15) << " " << XT::Common::to_string(std::abs(diff), 15) << std::endl;
    entropy_file.close();
    last_entropy_ = entropy;
  }

  template <class GV, class LocalDiscreteFunctionType, class MergedQuadratureType, class BasisValsType>
  double compute_entropy(LocalDiscreteFunctionType& local_u,
                         const GV& grid_view,
                         const MergedQuadratureType& merged_quads,
                         const BasisValsType& basis_vals,
                         std::vector<double>& vals)
  {
    vals.clear();
    XT::Common::FieldVector<RangeFieldType, dimRange> local_u_vec;
    for (auto&& entity : Dune::elements(grid_view)) {
      local_u->bind(entity);
      for (size_t jj = 0; jj < dimRange; ++jj)
        local_u_vec[jj] = local_u->dofs().get_entry(jj);
      const auto alpha = op_.entropy_solver().entropy_flux().get_alpha(local_u_vec, true)->first;
      size_t kk = 0;
      for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it, ++kk) {
        const auto& quad_point = *it;
        const auto alpha_n_b = alpha * basis_vals[kk];
        vals.push_back(std::exp(alpha_n_b) * (alpha_n_b - 1) * quad_point.weight());
      } // quad_points
    } // entities
    return compensated_sum(vals);
  }


  double compensated_sum(const std::vector<double>& input) const
  {
    double sum = 0.0;
    double c = 0.0; // A running compensation for lost low-order bits.

    for (const double val : input) {
      double t = sum + val;
      if (std::abs(sum) >= std::abs(val))
        c += (sum - t) + val; // If sum is bigger, low-order digits of input[i] are lost.
      else
        c += (val - t) + sum; // Else low-order digits of sum are lost.
      sum = t;
    }
    return sum + c; // Correction only applied once in the very end.
  }


  const OperatorType& op_;
  const RangeFieldType r_;
  std::unique_ptr<DiscreteFunctionType> u_i_;
  const MatrixType A_;
  const VectorType b_;
  const VectorType c_;
  std::vector<std::unique_ptr<DiscreteFunctionType>> stages_k_;
  const size_t num_stages_;
  std::vector<double> relaxationupdate_stages_;
  double last_entropy_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_EXPLICIT_RUNGEKUTTA_HH
