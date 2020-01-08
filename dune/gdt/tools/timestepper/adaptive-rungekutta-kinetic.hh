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

#ifndef DUNE_GDT_TIMESTEPPER_KINETIC_ADAPTIVE_RUNGEKUTTA_HH
#define DUNE_GDT_TIMESTEPPER_KINETIC_ADAPTIVE_RUNGEKUTTA_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/string.hh>

#include "adaptive-rungekutta.hh"


namespace Dune {
namespace GDT {


/** \brief Time stepper using adaptive Runge Kutta methods
 *
 * Timestepper using adaptive Runge Kutta methods to solve equations of the form u_t = r * L(u, t) where u is a
 * discrete function, L an operator acting on u and \alpha a scalar factor (e.g. -1).
 * The specific Runge Kutta method can be chosen as the third template argument. If your desired Runge Kutta method is
 * not contained in AdaptiveRungeKuttaMethods, choose AdaptiveRungeKuttaMethods::other and
 * supply a DynamicMatrix< RangeFieldType > A and vectors b_1, b_2 (DynamicVector< RangeFieldType >) and c
 * (DynamicVector< RangeFieldType >) in the constructor. Here, A, b_1, b_2 and c form the butcher tableau (see
 * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Embedded_methods, A is composed of the coefficients
 * a_{ij}, b_1 of b_j, b_2 of b_j^* and c of c_j). The default is the Dormand-Prince RK45 method.
 * In each time step, the error is estimated using the difference between the two solutions obtained using either b_1 or
 * b_2. If the estimated error is higher than a specified tolerance tol, the calculation is repeated with a smaller
 * time step. The tolerance tol and the error estimate are also used to estimate the optimal time step length for the
 * next time step via dt_new = dt_old*min(max(0.9*(tol/error)^(1/5), scale_factor_min), scale_factor_max_);
 *
 * Notation: For an s-stage method,
 * \mathbf{u}^{n+1} = \mathbf{u}^n + dt \sum_{i=0}^{s-1} b_i \mathbf{k}_i
 * \mathbf{k}_i = L(\mathbf{u}_i, t^n + dt c_i)
 * \mathbf{u}_i = \mathbf{u}^n + dt \sum_{j=0}^{i-1} a_{ij} \mathbf{k}_j,
 *
 * \tparam OperatorImp Type of operator L
 * \tparam DiscreteFunctionImp Type of initial values and solution at a fixed time
 * \tparam method Adaptive Runge-Kutta method that is used (default is AdaptiveRungeKuttaMethods::dormand_prince)
 */
template <class OperatorImp,
          class DiscreteFunctionImp,
          class EntropyFluxImp,
          TimeStepperMethods method = TimeStepperMethods::dormand_prince>
class KineticAdaptiveRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  typedef TimeStepperInterface<DiscreteFunctionImp> BaseType;
  typedef typename internal::AdaptiveButcherArrayProvider<typename BaseType::RangeFieldType, method>
      ButcherArrayProviderType;

public:
  typedef OperatorImp OperatorType;
  typedef DiscreteFunctionImp DiscreteFunctionType;

  typedef typename DiscreteFunctionType::DomainFieldType DomainFieldType;
  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename Dune::DynamicMatrix<RangeFieldType> MatrixType;
  typedef typename Dune::DynamicVector<RangeFieldType> VectorType;
  typedef typename std::vector<std::pair<RangeFieldType, DiscreteFunctionType>> SolutionType;
  static constexpr size_t q = ButcherArrayProviderType::q;
  using EntropyFluxType = EntropyFluxImp;
  using BaseType::dimRange;

  DiscreteFunctionType get_u(const DiscreteFunctionType& alpha)
  {
    DiscreteFunctionType u(alpha.space());
    const auto local_u = u.local_discrete_function();
    const auto local_alpha = alpha.local_discrete_function();
    XT::Common::FieldVector<RangeFieldType, dimRange> local_alpha_vec;
    XT::Common::FieldVector<RangeFieldType, dimRange> local_u_vec;
    for (auto&& entity : Dune::elements(alpha.space().grid_view())) {
      local_alpha->bind(entity);
      local_u->bind(entity);
      for (size_t ii = 0; ii < dimRange; ++ii)
        local_alpha_vec[ii] = local_alpha->dofs().get_entry(ii);
      local_u_vec = entropy_flux_.get_u(local_alpha_vec);
      for (size_t ii = 0; ii < dimRange; ++ii)
        local_u->dofs().set_entry(ii, local_u_vec[ii]);
    }
    return u;
  }

  /**
   * \brief Constructor for AdaptiveRungeKuttaTimeStepper time stepper
   * \param op Operator L
   * \param initial_values Discrete function containing initial values for u at time t_0.
   * \param r Scalar factor (see above, default is 1)
   * \param t_0 Initial time (default is 0)
   * \param tol Error tolerance for the adaptive scheme (default is 1e-4)
   * \param scale_factor_min Minimum allowed factor for time step scaling (default is 0.2).
   * \param scale_factor_max Maximum allowed factor for time step scaling (default is 5).
   * \param A Coefficient matrix (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param b_1 First set of coefficients (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param b_2 Second set of coefficients (only provide if you use AdaptiveRungeKuttaMethods::other)
   * \param c Coefficients for time steps (only provide if you use AdaptiveRungeKuttaMethods::other)
   */
  KineticAdaptiveRungeKuttaTimeStepper(const OperatorType& op,
                                       const EntropyFluxType& entropy_flux,
                                       DiscreteFunctionType& initial_values,
                                       const RangeFieldType r = 1.0,
                                       const double t_0 = 0.0,
                                       const RangeFieldType tol = 1e-3,
                                       const RangeFieldType scale_factor_min = 0.2,
                                       const RangeFieldType scale_factor_max = 5,
                                       const MatrixType& A = ButcherArrayProviderType::A(),
                                       const VectorType& b_1 = ButcherArrayProviderType::b_1(),
                                       const VectorType& b_2 = ButcherArrayProviderType::b_2(),
                                       const VectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , op_(op)
    , entropy_flux_(entropy_flux)
    , r_(r)
    , atol_(tol)
    , rtol_(tol)
    , scale_factor_min_(scale_factor_min)
    , scale_factor_max_(scale_factor_max)
    , alpha_tmp_(BaseType::current_solution())
    , alpha_np1_(BaseType::current_solution())
    , A_(A)
    , b_1_(b_1)
    , b_2_(b_2)
    , c_(c)
    , b_diff_(b_2_ - b_1_)
    , num_stages_(A_.rows())
  {
    assert(Dune::XT::Common::FloatCmp::ge(atol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::ge(rtol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::le(scale_factor_min_, 1.0));
    assert(Dune::XT::Common::FloatCmp::ge(scale_factor_max_, 1.0));
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_1_.size() == A_.rows());
    assert(b_2_.size() == A_.rows());
    assert(c_.size() == A_.rows());
#ifndef NDEBUG
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        assert(Dune::XT::Common::FloatCmp::eq(A_[ii][jj], 0.0)
               && "A has to be a lower triangular matrix with 0 on the main diagonal");
      }
    }
#endif // NDEBUG
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(current_solution());
    }
  } // constructor AdaptiveRungeKuttaTimeStepper

  using BaseType::current_solution;
  using BaseType::current_time;
  using BaseType::solve;

  virtual RangeFieldType solve(const RangeFieldType t_end,
                               const RangeFieldType initial_dt,
                               const size_t num_save_steps,
                               const size_t num_output_steps,
                               const bool save_solution,
                               const bool visualize,
                               const bool write_discrete,
                               const bool write_exact,
                               const std::string prefix,
                               typename BaseType::DiscreteSolutionType& sol,
                               const typename BaseType::VisualizerType& visualizer,
                               const typename BaseType::StringifierType& stringifier,
                               const typename BaseType::GridFunctionType& exact_solution) override final
  {
    const auto ret = BaseType::solve(t_end,
                                     initial_dt,
                                     num_save_steps,
                                     num_output_steps,
                                     save_solution,
                                     visualize,
                                     write_discrete,
                                     write_exact,
                                     prefix,
                                     sol,
                                     visualizer,
                                     stringifier,
                                     exact_solution);
    // in a fractional step scheme, we cannot use last_stage_of_previous_step
    last_stage_of_previous_step_ = nullptr;
    return ret;
  }

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    RangeFieldType actual_dt = std::min(dt, max_dt);
    RangeFieldType mixed_error = std::numeric_limits<RangeFieldType>::max();
    RangeFieldType time_step_scale_factor = 1.0;
    const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1};
    auto r_it = r_sequence.begin();

    auto& t = current_time();
    auto& alpha_n = current_solution();

    while (mixed_error > 1.) {
      bool skip_error_computation = false;
      actual_dt *= time_step_scale_factor;
      size_t first_stage_to_compute = 0;
      if (last_stage_of_previous_step_) {
        stages_k_[0].dofs().vector() = last_stage_of_previous_step_->dofs().vector();
        first_stage_to_compute = 1;
      }

      bool consider_regularization = actual_dt < 1e-5 ? true : false;
      for (size_t ii = first_stage_to_compute; ii < num_stages_; ++ii) {
        stages_k_[ii].dofs().vector() *= 0.;
        alpha_tmp_.dofs().vector() = alpha_n.dofs().vector();
        for (size_t jj = 0; jj < ii; ++jj)
          alpha_tmp_.dofs().vector() += stages_k_[jj].dofs().vector() * (actual_dt * r_ * (A_[ii][jj]));
        try {
          op_.apply(alpha_tmp_.dofs().vector(),
                    stages_k_[ii].dofs().vector(),
                    XT::Common::Parameter(
                        {{"t", {t + actual_dt * c_[ii]}}, {"reg", {static_cast<double>(consider_regularization)}}}));
          const auto& reg_indicators = op_.reg_indicators();
          if (consider_regularization
              && std::any_of(reg_indicators.begin(), reg_indicators.end(), [](const bool& a) { return a; })) {
            const auto& grid_view = alpha_n.space().grid_view();
            const auto local_alpha = alpha_n.local_discrete_function();
            const auto& basis_functions = entropy_flux_.basis_functions();
            const auto alpha_one = basis_functions.alpha_one();
            XT::Common::FieldVector<RangeFieldType, dimRange> local_alpha_vec;
            if (++r_it == r_sequence.end())
              DUNE_THROW(Dune::InvalidStateException, "Fully regularized, still fails!");
            const auto r = *r_it;
            for (auto&& entity : Dune::elements(grid_view)) {
              const auto entity_index = grid_view.indexSet().index(entity);
              if (reg_indicators[entity_index]) {
                std::cout << "Regularized on entity " << entity_index << " with r = " << r << std::endl;
                local_alpha->bind(entity);
                for (size_t jj = 0; jj < dimRange; ++jj)
                  local_alpha_vec[jj] = local_alpha->dofs().get_entry(jj);
                const auto old_density = basis_functions.density(entropy_flux_.get_u(local_alpha_vec));
                const auto alpha_iso = basis_functions.alpha_iso(old_density);
                for (size_t jj = 0; jj < dimRange; ++jj)
                  local_alpha_vec[jj] = (1 - r) * local_alpha_vec[jj] + r * alpha_iso[jj];
                const auto reg_density = basis_functions.density(entropy_flux_.get_u(local_alpha_vec));
                const auto factor = std::log(old_density / reg_density);
                for (size_t jj = 0; jj < dimRange; ++jj) {
                  local_alpha_vec[jj] += alpha_one[jj] * factor;
                  local_alpha->dofs().set_entry(jj, local_alpha_vec[jj]);
                }
              }
            }
            last_stage_of_previous_step_ = nullptr;
            DUNE_THROW(Dune::MathError, "Regularization done!");
          }
        } catch (const Dune::MathError& e) {
          mixed_error = 1e10;
          skip_error_computation = true;
          time_step_scale_factor = 0.5;
          // std::cout << "MathError! " << e.what() << std::endl;
          break;
#if HAVE_TBB
        } catch (const tbb::captured_exception& e) {
          mixed_error = 1e10;
          skip_error_computation = true;
          time_step_scale_factor = 0.5;
          // std::cout << "TBB error! " << e.what() << std::endl;
          break;
#endif
        }
      }

      if (!skip_error_computation) {
        // calculate two approximations of alpha at timestep n+1.
        alpha_tmp_.dofs().vector() = alpha_n.dofs().vector();
        alpha_np1_.dofs().vector() = alpha_n.dofs().vector();
        for (size_t ii = 0; ii < num_stages_; ++ii) {
          alpha_np1_.dofs().vector().axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii].dofs().vector());
          alpha_tmp_.dofs().vector().axpy(actual_dt * r_ * b_2_[ii], stages_k_[ii].dofs().vector());
        }

        // calculate error
        const auto* alpha_tmp_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_tmp_.dofs().vector());
        const auto* alpha_np1_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_np1_.dofs().vector());
        mixed_error =
            XT::Common::transform_reduce(alpha_tmp_data,
                                         alpha_tmp_data + alpha_tmp_.dofs().vector().size(),
                                         alpha_np1_data,
                                         0.,
                                         /*reduction*/ [](const auto& a, const auto& b) { return std::max(a, b); },
                                         /*transformation*/
                                         [atol = atol_, rtol = rtol_](const auto& a, const auto& b) {
                                           return std::abs(a - b) / (atol + std::max(std::abs(a), std::abs(b)) * rtol);
                                         });
        // scale dt to get the estimated optimal time step length
        time_step_scale_factor =
            std::min(std::max(0.8 * std::pow(1. / mixed_error, 1. / (q + 1.)), scale_factor_min_), scale_factor_max_);
      }
    } // while (mixed_error > 1.)

    alpha_n.dofs().vector() = alpha_np1_.dofs().vector();

    // if (!last_stage_of_previous_step_)
    // last_stage_of_previous_step_ = Dune::XT::Common::make_unique<DiscreteFunctionType>(alpha_n);
    // last_stage_of_previous_step_->dofs().vector() = stages_k_[num_stages_ - 1].dofs().vector();

    t += actual_dt;

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

private:
  const OperatorType& op_;
  const EntropyFluxType& entropy_flux_;
  const RangeFieldType r_;
  const RangeFieldType atol_;
  const RangeFieldType rtol_;
  const RangeFieldType scale_factor_min_;
  const RangeFieldType scale_factor_max_;
  DiscreteFunctionType alpha_tmp_;
  DiscreteFunctionType alpha_np1_;
  const MatrixType A_;
  const VectorType b_1_;
  const VectorType b_2_;
  const VectorType c_;
  const VectorType b_diff_;
  std::vector<DiscreteFunctionType> stages_k_;
  const size_t num_stages_;
  std::unique_ptr<DiscreteFunctionType> last_stage_of_previous_step_;
}; // class KineticAdaptiveRungeKuttaTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_KINETIC_ADAPTIVE_RUNGEKUTTA_HH
