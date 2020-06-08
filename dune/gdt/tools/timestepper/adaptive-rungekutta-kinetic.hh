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
#include <dune/xt/common/numeric.hh>
#include <dune/xt/common/string.hh>

#include "adaptive-rungekutta.hh"


namespace Dune {
namespace GDT {


/** \brief Specialization of AdaptiveRungeKuttaTimestepper for coordinate-transformed minimum entropy moment models
 *
 * Specializations include:
 *  - catches exceptions that are thrown and tries again with reduced timestep instead of aborting
 *  - a step to ensure a minimum density or apply other regularization
 *  - isotropic regularization (unfinished and not really tested, do not use unless you know what you are doing)
 */
template <class OperatorImp,
          class MinDensitySetterType,
          class DiscreteFunctionImp,
          class EntropyFluxImp,
          TimeStepperMethods method = TimeStepperMethods::dormand_prince>
class KineticAdaptiveRungeKuttaTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  using BaseType = TimeStepperInterface<DiscreteFunctionImp>;
  using ButcherArrayProviderType =
      typename internal::AdaptiveButcherArrayProvider<typename BaseType::RangeFieldType, method>;

public:
  using OperatorType = OperatorImp;
  using DiscreteFunctionType = DiscreteFunctionImp;
  using typename BaseType::DiscreteSolutionType;
  using DomainFieldType = typename DiscreteFunctionType::DomainFieldType;
  using RangeFieldType = typename DiscreteFunctionType::RangeFieldType;
  using MatrixType = typename Dune::DynamicMatrix<RangeFieldType>;
  using VectorType = typename Dune::DynamicVector<RangeFieldType>;
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

  KineticAdaptiveRungeKuttaTimeStepper(const OperatorType& op,
                                       const MinDensitySetterType& min_density_setter,
                                       const EntropyFluxType& entropy_flux,
                                       DiscreteFunctionType& initial_values,
                                       const bool /*use_first_same_as_last_property*/ = true,
                                       const RangeFieldType r = 1.0,
                                       const RangeFieldType atol = 1e-3,
                                       const RangeFieldType rtol = 1e-2,
                                       const double t_0 = 0.0,
                                       const RangeFieldType scale_factor_min = 0.2,
                                       const RangeFieldType scale_factor_max = 5,
                                       const MatrixType& A = ButcherArrayProviderType::A(),
                                       const VectorType& b_1 = ButcherArrayProviderType::b_1(),
                                       const VectorType& b_2 = ButcherArrayProviderType::b_2(),
                                       const VectorType& c = ButcherArrayProviderType::c())
    : BaseType(t_0, initial_values)
    , min_density_setter_(min_density_setter)
    , op_(op)
    , entropy_flux_(entropy_flux)
    , r_(r)
    , atol_(atol)
    , rtol_(rtol)
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
    , initial_values_evaluated_(false)
  {
    assert(Dune::XT::Common::FloatCmp::ge(atol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::ge(rtol_, 0.0));
    assert(Dune::XT::Common::FloatCmp::le(scale_factor_min_, 1.0));
    assert(Dune::XT::Common::FloatCmp::ge(scale_factor_max_, 1.0));
    assert(A_.rows() == A_.cols() && "A has to be a square matrix");
    assert(b_1_.size() == A_.rows());
    assert(b_2_.size() == A_.rows());
    assert(c_.size() == A_.rows());
    for (size_t ii = 0; ii < A_.rows(); ++ii) {
      for (size_t jj = ii; jj < A_.cols(); ++jj) {
        DUNE_THROW_IF(XT::Common::FloatCmp::ne(A_[ii][jj], 0.0),
                      XT::Common::Exceptions::wrong_input_given,
                      "A has to be a lower triangular matrix with 0 on the main diagonal");
      }
    }
    // store as many discrete functions as needed for intermediate stages
    for (size_t ii = 0; ii < num_stages_; ++ii) {
      stages_k_.emplace_back(current_solution());
    }
  } // constructor KineticAdaptiveRungeKuttaTimeStepper

  using BaseType::current_solution;
  using BaseType::current_time;
  using BaseType::solve;

  bool regularize_if_needed(const bool /*consider_regularization*/,
                            typename std::vector<RangeFieldType>::const_iterator& /*r_it*/,
                            const std::vector<RangeFieldType>& /*r_sequence*/)
  {
    return false;
#if 0
    const auto& reg_indicators = op_.reg_indicators();
    auto& alpha_n = current_solution();
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
          local_alpha->bind(entity);
          for (size_t jj = 0; jj < dimRange; ++jj)
            local_alpha_vec[jj] = local_alpha->dofs().get_entry(jj);
          std::cout << "Regularized on entity " << entity_index << " with r = " << r << "and alpha "
                    << XT::Common::to_string(local_alpha_vec) << std::endl;
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
      return true;
    }
    return false;
#endif
  } // void regularize_if_needed(...)

  void set_op_param(const std::string& key, const RangeFieldType& val)
  {
    static std::vector<RangeFieldType> tmp_vec(1);
    tmp_vec[0] = val;
    op_param_.set(key, tmp_vec, true);
  }

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    RangeFieldType actual_dt = std::min(dt, max_dt);
    // regularization is currently unused but may be used in the near future
    static const bool consider_regularization = false;
    // set_op_param("dt", actual_dt);
    RangeFieldType mixed_error = std::numeric_limits<RangeFieldType>::max();
    RangeFieldType time_step_scale_factor = 1.0;
    static const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1};
    auto r_it = r_sequence.begin();

    auto& t = current_time();
    auto& alpha_n = current_solution();
    const auto num_dofs = alpha_n.dofs().vector().size();
    size_t first_stage_to_compute = 0;
    if (first_same_as_last_ && last_stage_of_previous_step_) {
      stages_k_[0].dofs().vector() = last_stage_of_previous_step_->dofs().vector();
      first_stage_to_compute = 1;
    }
    first_same_as_last_ = true;
    while (!(mixed_error < 1.)) {
      bool skip_error_computation = false;
      actual_dt *= time_step_scale_factor;
      for (size_t ii = first_stage_to_compute; ii < num_stages_ - 1; ++ii) {
        // set_op_param("t", t + actual_dt * c_[ii]);
        std::fill_n(&(stages_k_[ii].dofs().vector()[0]), num_dofs, 0.);
        alpha_tmp_.dofs().vector() = alpha_n.dofs().vector();
        for (size_t jj = 0; jj < ii; ++jj)
          alpha_tmp_.dofs().vector().axpy(actual_dt * r_ * A_[ii][jj], stages_k_[jj].dofs().vector());
        try {
          op_.apply(alpha_tmp_.dofs().vector(), stages_k_[ii].dofs().vector(), op_param_);
          if (regularize_if_needed(consider_regularization, r_it, r_sequence)) {
            mixed_error = 1e10;
            skip_error_computation = true;
            time_step_scale_factor = 0.9;
            break;
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
      } // stages

      if (!skip_error_computation) {
        // compute alpha^{n+1}
        alpha_np1_.dofs().vector() = alpha_n.dofs().vector();
        for (size_t ii = 0; ii < num_stages_ - 1; ++ii)
          alpha_np1_.dofs().vector().axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii].dofs().vector());

        // calculate last stage
        // set_op_param("t", t + actual_dt * c_[num_stages_ - 1]);
        std::fill_n(&(stages_k_[num_stages_ - 1].dofs().vector()[0]), num_dofs, 0.);
        try {
          op_.apply(alpha_np1_.dofs().vector(), stages_k_[num_stages_ - 1].dofs().vector(), op_param_);
          if (regularize_if_needed(consider_regularization, r_it, r_sequence)) {
            mixed_error = 1e10;
            time_step_scale_factor = 0.9;
            continue;
          }
        } catch (const Dune::MathError& e) {
          mixed_error = 1e10;
          time_step_scale_factor = 0.5;
          // std::cout << "MathError! " << e.what() << std::endl;
          continue;
#if HAVE_TBB
        } catch (const tbb::captured_exception& e) {
          mixed_error = 1e10;
          time_step_scale_factor = 0.5;
          continue;
#endif
        }

        // calculate second approximations of alpha at timestep n+1.
        alpha_tmp_.dofs().vector() = alpha_n.dofs().vector();
        for (size_t ii = 0; ii < num_stages_; ++ii)
          alpha_tmp_.dofs().vector().axpy(actual_dt * r_ * b_2_[ii], stages_k_[ii].dofs().vector());

        const auto* alpha_tmp_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_tmp_.dofs().vector());
        const auto* alpha_np1_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_np1_.dofs().vector());

        // if b == NaN, std::max(a, b) might return a, so mixed_error might be non-NaN though there are NaNs in the
        // vectors So check for NaNs before
        bool nan_found = check_for_nan(alpha_tmp_data, alpha_np1_data, num_dofs);
        if (nan_found) {
          mixed_error = 1e10;
          time_step_scale_factor = 0.5;
        } else {
          // calculate error
          mixed_error = XT::Common::transform_reduce(
              alpha_tmp_data,
              alpha_tmp_data + num_dofs,
              alpha_np1_data,
              0.,
              /*reduction*/ [](const auto& a, const auto& b) { return std::max(a, b); },
              /*transformation*/
              [atol = atol_, rtol = rtol_](const auto& a, const auto& b) {
                return std::abs(a - b) / (atol + std::max(std::abs(a), std::abs(b)) * rtol);
              });
          // std::cout << mixed_error << std::endl;
          // scale dt to get the estimated optimal time step length
          time_step_scale_factor =
              std::min(std::max(0.8 * std::pow(1. / mixed_error, 1. / (q + 1.)), scale_factor_min_), scale_factor_max_);

          // maybe adjust alpha to enforce a minimum density or avoid problems with matrix conditions
          if (mixed_error < 1.
              && min_density_setter_.apply_with_dt(alpha_np1_.dofs().vector(), alpha_np1_.dofs().vector(), actual_dt)) {
            // we cannot use the first-same-as-last property for the next time step if we changed alpha
            first_same_as_last_ = false;
          }
        }
      } // if (!skip_error_computation)
    } // while (mixed_error > 1.)

    alpha_n.dofs().vector() = alpha_np1_.dofs().vector();
    this->dts_.push_back(actual_dt);

    if (!last_stage_of_previous_step_)
      last_stage_of_previous_step_ = std::make_unique<DiscreteFunctionType>(alpha_n);
    last_stage_of_previous_step_->dofs().vector() = stages_k_[num_stages_ - 1].dofs().vector();

    if (operator_evaluations_) {
      for (size_t ii = first_stage_to_compute; ii < num_stages_; ++ii)
        operator_evaluations_->push_back(stages_k_[ii].dofs().vector());
    }

    t += actual_dt;

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

  RangeFieldType next_n_steps(const size_t n,
                              const RangeFieldType t_end,
                              const RangeFieldType initial_dt,
                              const bool /*output_progress*/,
                              const bool /*with_half_steps*/,
                              DiscreteSolutionType& sol) override final
  {
    RangeFieldType dt = initial_dt;
    RangeFieldType t = current_time();
    assert(XT::Common::FloatCmp::ge(t_end - t, 0.0));
    size_t time_step_counter = 0;

    if (!initial_values_evaluated_ && XT::Common::FloatCmp::eq(t, 0.0)) {
      sol.insert(sol.end(), std::make_pair(t, current_solution()));
      ++time_step_counter;
      initial_values_evaluated_ = true;
    }

    while (XT::Common::FloatCmp::lt(t, t_end) && time_step_counter < n) {
      RangeFieldType max_dt = XT::Common::FloatCmp::ge(t + dt, t_end) ? t_end - t : dt;
      dt = step(dt, max_dt);
      t = current_time();
      sol.insert(sol.end(), std::make_pair(t, current_solution()));
      ++time_step_counter;
    }
    return dt;
  } // ... next_n_steps(...)

private:
  bool check_for_nan(const RangeFieldType* vec1, const RangeFieldType* vec2, const size_t num_dofs)
  {
    for (size_t ii = 0; ii < num_dofs; ++ii)
      if (std::isnan(vec1[ii] + vec2[ii]))
        return true;
    return false;
  }

  const MinDensitySetterType min_density_setter_;
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
  bool first_same_as_last_;
  XT::Common::Parameter op_param_;
  bool initial_values_evaluated_;
  using BaseType::operator_evaluations_;
}; // class KineticAdaptiveRungeKuttaTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_KINETIC_ADAPTIVE_RUNGEKUTTA_HH
