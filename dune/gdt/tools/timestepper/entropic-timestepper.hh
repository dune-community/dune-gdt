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

#ifndef DUNE_GDT_TIMESTEPPER_ENTROPIC_HH
#define DUNE_GDT_TIMESTEPPER_ENTROPIC_HH

#include <utility>

#include <dune/gdt/operators/interfaces.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/common/numeric.hh>
#include <dune/xt/common/string.hh>

#include "adaptive-rungekutta.hh"

namespace Dune {
namespace GDT {


template <class OperatorImp,
          class EntropyFluxImp,
          class ValuesImp,
          TimeStepperMethods method = TimeStepperMethods::dormand_prince>
class EntropicTimeStepper
{
  using OperatorType = OperatorImp;
  using EntropyFluxType = EntropyFluxImp;
  using ValuesType = ValuesImp;
  using RangeFieldType = double;
  static constexpr size_t dimRange = ValuesType::dimRange;
  using VisualizerType = XT::Functions::VisualizerInterface<dimRange, 1, RangeFieldType>;
  using CoeffMatType = Dune::DynamicMatrix<RangeFieldType>;
  using CoeffVecType = Dune::DynamicVector<RangeFieldType>;
  using ButcherArrayProviderType = internal::AdaptiveButcherArrayProvider<RangeFieldType, method>;
  static constexpr size_t q = ButcherArrayProviderType::q;

public:
  EntropicTimeStepper(OperatorType& op,
                      const EntropyFluxType& entropy_flux,
                      ValuesType& alpha_initial,
                      const RangeFieldType r = 1.0,
                      const double t_0 = 0.0,
                      const RangeFieldType tol = 1e-3,
                      const RangeFieldType scale_factor_min = 0.2,
                      const RangeFieldType scale_factor_max = 5,
                      const CoeffMatType& A = ButcherArrayProviderType::A(),
                      const CoeffVecType& b_1 = ButcherArrayProviderType::b_1(),
                      const CoeffVecType& b_2 = ButcherArrayProviderType::b_2(),
                      const CoeffVecType& c = ButcherArrayProviderType::c())
    : t_(t_0)
    , op_(op)
    , entropy_flux_(entropy_flux)
    , alpha_(alpha_initial)
    , alpha_tmp_(alpha_)
    , alpha_np1_(alpha_)
    , r_(r)
    , atol_(tol)
    , rtol_(tol)
    , scale_factor_min_(scale_factor_min)
    , scale_factor_max_(scale_factor_max)
    , A_(A)
    , b_1_(b_1)
    , b_2_(b_2)
    , c_(c)
    , b_diff_(b_2_ - b_1_)
    , num_stages_(A_.rows())
  {
    assert(Dune::XT::Common::FloatCmp::gt(tol, 0.0));
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
      stages_k_.emplace_back(alpha_);
    }
  } // constructor AdaptiveRungeKuttaTimeStepper

  const RangeFieldType& current_time() const
  {
    return t_;
  }

  const ValuesType& current_solution() const
  {
    return alpha_;
  }

  RangeFieldType solve(const RangeFieldType t_end,
                       const RangeFieldType initial_dt,
                       const size_t num_save_steps,
                       const size_t num_output_steps,
                       const bool visualize,
                       const std::string prefix,
                       const VisualizerType& visualizer)
  {
    RangeFieldType dt = initial_dt;
    RangeFieldType t = current_time();
    assert(Dune::XT::Common::FloatCmp::ge(t_end, t));
    size_t time_step_counter = 0;

    const RangeFieldType save_interval = (t_end - t) / num_save_steps;
    const RangeFieldType output_interval =
        (num_output_steps == 0 ? std::numeric_limits<RangeFieldType>::max() : (t_end - t) / num_output_steps);
    RangeFieldType next_save_time = t + save_interval > t_end ? t_end : t + save_interval;
    RangeFieldType next_output_time = t + output_interval > t_end ? t_end : t + output_interval;
    size_t save_step_counter = 1;

    // save/visualize initial solution
    write_files(visualize, current_solution(), prefix, 0, visualizer);

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
        write_files(visualize, current_solution(), prefix, save_step_counter, visualizer);
        next_save_time += save_interval;
        ++save_step_counter;
      }
      if (num_output_steps && (Dune::XT::Common::FloatCmp::ge(t, next_output_time) || num_output_steps == size_t(-1))) {
        std::cout << "time step " << time_step_counter << " done, time =" << t << ", current dt= " << dt << std::endl;
        next_output_time += output_interval;
      }
    } // while (t < t_end)

    return dt;
  } // ... solve(...)

  static void write_files(const bool visualize,
                          const ValuesType& alpha,
                          const std::string& prefix,
                          const size_t step,
                          const VisualizerType& visualizer)
  {
    if (visualize)
      alpha.visualize(prefix + "_" + XT::Common::to_string(step), visualizer);
  }

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt)
  {
    RangeFieldType actual_dt = std::min(dt, max_dt);
    RangeFieldType mixed_error = std::numeric_limits<RangeFieldType>::max();
    RangeFieldType time_step_scale_factor = 1.0;
    const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1};
    auto r_it = r_sequence.begin();

    auto& t = t_;
    auto& alpha_n = alpha_;

    size_t count = 0;
    while (mixed_error > 1.) {
      if (count++ != 0)
        std::cout << "mixed_error = " << mixed_error << std::endl;
      bool skip_error_computation = false;
      actual_dt *= time_step_scale_factor;
      size_t first_stage_to_compute = 0;
      if (last_stage_of_previous_step_) {
        stages_k_[0].values() = last_stage_of_previous_step_->values();
        first_stage_to_compute = 1;
      }

      bool consider_regularization = actual_dt < 1e-5 ? true : false;
      for (size_t ii = first_stage_to_compute; ii < num_stages_; ++ii) {
        stages_k_[ii].set_zero();
        alpha_tmp_.copy_values(alpha_n);
        for (size_t jj = 0; jj < ii; ++jj)
          alpha_tmp_.axpy(actual_dt * r_ * (A_[ii][jj]), stages_k_[jj]);
        try {
          op_.apply(alpha_tmp_, stages_k_[ii], consider_regularization);
          const auto& reg_indicators = op_.reg_indicators();
          if (consider_regularization
              && std::any_of(reg_indicators.begin(), reg_indicators.end(), [](const bool& a) { return a; })) {
            std::cout << "Regularization considered" << std::endl;
            if (++r_it == r_sequence.end())
              DUNE_THROW(Dune::InvalidStateException, "Fully regularized, still fails!");
            const auto r = *r_it;
            op_.regularize(r, alpha_n);
            last_stage_of_previous_step_ = nullptr;
            DUNE_THROW(Dune::MathError, "Regularization done!");
          }
        } catch (const Dune::MathError& e) {
          std::cout << "Math error caught" << e.what() << std::endl;
          mixed_error = 1e10;
          skip_error_computation = true;
          time_step_scale_factor = 0.1;
          break;
#if HAVE_TBB
        } catch (const tbb::captured_exception& e) {
          std::cout << "TBB exception caught" << std::endl;
          mixed_error = 1e10;
          skip_error_computation = true;
          time_step_scale_factor = 0.1;
          break;
#endif
        }
      }

      if (!skip_error_computation) {
        // calculate two approximations of alpha at timestep n+1.
        alpha_tmp_.copy_values(alpha_n);
        alpha_np1_.copy_values(alpha_n);
        for (size_t ii = 0; ii < num_stages_; ++ii) {
          alpha_np1_.axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii]);
          alpha_tmp_.axpy(actual_dt * r_ * b_2_[ii], stages_k_[ii]);
        }
        // ensure that we do not fall below vacuum density
        const auto& basis_functions = op_.entropy_flux().basis_functions();
        const auto psi_min = op_.psi_min();
        for (size_t nn = 0; nn < alpha_np1_.size(); ++nn) {
          basis_functions.adjust_alpha_to_ensure_min_density(alpha_tmp_[nn], psi_min);
          basis_functions.adjust_alpha_to_ensure_min_density(alpha_np1_[nn], psi_min);
        }
        // calculate error
        mixed_error =
            XT::Common::transform_reduce(alpha_tmp_.data(),
                                         alpha_tmp_.data() + alpha_tmp_.num_elements(),
                                         alpha_np1_.data(),
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
    } // while (mixed_error > tol_)

    alpha_n.copy_values(alpha_np1_);

    t += actual_dt;

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

private:
  double t_;
  OperatorType& op_;
  const EntropyFluxType& entropy_flux_;
  ValuesType& alpha_;
  ValuesType alpha_tmp_;
  ValuesType alpha_np1_;
  const RangeFieldType r_;
  const RangeFieldType atol_;
  const RangeFieldType rtol_;
  const RangeFieldType scale_factor_min_;
  const RangeFieldType scale_factor_max_;
  const CoeffMatType A_;
  const CoeffVecType b_1_;
  const CoeffVecType b_2_;
  const CoeffVecType c_;
  const CoeffVecType b_diff_;
  std::vector<ValuesType> stages_k_;
  const size_t num_stages_;
  std::unique_ptr<ValuesType> last_stage_of_previous_step_;
}; // class KineticAdaptiveRungeKuttaTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_ENTROPIC_HH
