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

#if HAVE_TBB && __has_include(<tbb/tbb_exception.h>)
#  include <tbb/tbb_exception.h>
#endif

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

  using DomainFieldType = typename DiscreteFunctionType::DomainFieldType;
  using RangeFieldType = typename DiscreteFunctionType::RangeFieldType;
  using MatrixType = typename Dune::DynamicMatrix<RangeFieldType>;
  using VectorType = typename Dune::DynamicVector<RangeFieldType>;
  using SolutionType = typename std::vector<std::pair<RangeFieldType, DiscreteFunctionType>>;
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
                                       const bool use_first_same_as_last_property = true,
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
    , alpha_tmp_(BaseType::current_solution().copy_as_discrete_function())
    , alpha_np1_(BaseType::current_solution().copy_as_discrete_function())
    , A_(A)
    , b_1_(b_1)
    , b_2_(b_2)
    , c_(c)
    , b_diff_(b_2_ - b_1_)
    , num_stages_(A_.rows())
    , first_same_as_last_(use_first_same_as_last_property)
    , relaxationupdate_stages_(num_stages_)
    , last_gamma_(1.)
    , last_entropy_(0)
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
      stages_k_.emplace_back(current_solution().copy_as_discrete_function());
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

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt, const std::string prefix) override final
  {
    RangeFieldType actual_dt = std::min(dt, max_dt);
    // regularization is currently unused but may be used in the near future
    static constexpr bool consider_regularization = false;
    set_op_param("dt", actual_dt);
    RangeFieldType mixed_error = std::numeric_limits<RangeFieldType>::max();
    RangeFieldType time_step_scale_factor = 1.0;
    static const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1};
    auto r_it = r_sequence.begin();
    auto& t = current_time();
    auto& alpha_n = current_solution();

    const auto local_alpha = alpha_n.local_discrete_function();
    std::vector<double> val_vector;
    const auto& grid_view = alpha_n.space().grid_view();
    const auto& quadratures = entropy_flux_.basis_functions().quadratures();
    const auto merged_quads = XT::Data::merged_quadrature(quadratures);
    const auto basis_vals = get_basis_vals(merged_quads);
    if (XT::Common::is_zero(t)) {
      last_entropy_ = compute_entropy(local_alpha, grid_view, merged_quads, basis_vals, val_vector);
      write_entropy(local_alpha, grid_view, 0, merged_quads, basis_vals, t, actual_dt, val_vector, prefix);
    }

    const auto num_dofs = alpha_n.dofs().vector().size();
    size_t first_stage_to_compute = 0;
    if (first_same_as_last_ && last_stage_of_previous_step_) {
      stages_k_[0]->dofs().vector() = last_stage_of_previous_step_->dofs().vector();
      relaxationupdate_stages_[0] = last_relaxation_of_previous_step_;
      first_stage_to_compute = 1;
    }
    first_same_as_last_ = true;
    const bool apply_gamma_relaxation = DXTC_CONFIG_GET("apply_gamma_relaxation", 0);
    // double min_gamma = 0.;
    double gamma = last_gamma_;
    double relaxationupdate(0.);
    while (!(mixed_error < 1.)) {
      relaxationupdate = 0.;
      gamma = last_gamma_;
      bool skip_error_computation = false;
      actual_dt *= time_step_scale_factor;
      set_op_param("dt", actual_dt);
      for (size_t ii = 0; ii < first_stage_to_compute; ++ii)
        relaxationupdate += b_1_[ii] * relaxationupdate_stages_[ii];
      for (size_t ii = first_stage_to_compute; ii < num_stages_ - 1; ++ii) {
        set_op_param("t", t + actual_dt * c_[ii]);
        std::fill_n(&(stages_k_[ii]->dofs().vector()[0]), num_dofs, 0.);
        alpha_tmp_->dofs().vector() = alpha_n.dofs().vector();
        for (size_t jj = 0; jj < ii; ++jj)
          alpha_tmp_->dofs().vector().axpy(actual_dt * r_ * A_[ii][jj], stages_k_[jj]->dofs().vector());
        try {
          relaxationupdate_stages_[ii] = 0.;
          op_.apply_and_compute_relax(
              alpha_tmp_->dofs().vector(), stages_k_[ii]->dofs().vector(), op_param_, relaxationupdate_stages_[ii]);
          relaxationupdate += b_1_[ii] * relaxationupdate_stages_[ii];
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
#if HAVE_TBB && __has_include(<tbb/tbb_exception.h>)
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
        alpha_np1_->dofs().vector() = alpha_n.dofs().vector();
        for (size_t ii = 0; ii < num_stages_ - 1; ++ii)
          alpha_np1_->dofs().vector().axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii]->dofs().vector());

        // calculate last stage
        set_op_param("t", t + actual_dt * c_[num_stages_ - 1]);
        std::fill_n(&(stages_k_[num_stages_ - 1]->dofs().vector()[0]), num_dofs, 0.);
        try {
          relaxationupdate_stages_[num_stages_ - 1] = 0.;
          op_.apply_and_compute_relax(alpha_np1_->dofs().vector(),
                                      stages_k_[num_stages_ - 1]->dofs().vector(),
                                      op_param_,
                                      relaxationupdate_stages_[num_stages_ - 1]);
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
#if HAVE_TBB && __has_include(<tbb/tbb_exception.h>)
        } catch (const tbb::captured_exception& e) {
          mixed_error = 1e10;
          time_step_scale_factor = 0.5;
          continue;
#endif
        }

        // calculate second approximations of alpha at timestep n+1.
        alpha_tmp_->dofs().vector() = alpha_n.dofs().vector();
        for (size_t ii = 0; ii < num_stages_; ++ii)
          alpha_tmp_->dofs().vector().axpy(actual_dt * r_ * b_2_[ii], stages_k_[ii]->dofs().vector());

        const auto* alpha_tmp_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_tmp_->dofs().vector());
        const auto* alpha_np1_data =
            XT::Common::VectorAbstraction<typename DiscreteFunctionType::VectorType>::data(alpha_np1_->dofs().vector());

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
              && min_density_setter_.apply_with_dt(
                  alpha_np1_->dofs().vector(), alpha_np1_->dofs().vector(), actual_dt)) {
            // we cannot use the first-same-as-last property for the next time step if we changed alpha
            first_same_as_last_ = false;
          }
        }
        if (mixed_error < 1.) {
          // visualize r(gamma)
          auto d = alpha_n.copy_as_discrete_function();
          d->dofs().vector() *= 0.;
          for (size_t ii = 0; ii < num_stages_ - 1; ++ii)
            d->dofs().vector().axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii]->dofs().vector());
          const auto local_d = d->local_discrete_function();
          // eta \circ eta_{\ast}^{\prime} = exp(x)(x - 1) for MaxwellBoltzmann entropy
          const bool visualize_gamma = DXTC_CONFIG_GET("visualize_gamma", 0);
          // std::cout << "relupdate: " << XT::Common::to_string(relaxationupdate, 15) << std::endl;
          if (visualize_gamma) {
            write_gamma(
                local_alpha, local_d, grid_view, merged_quads, basis_vals, actual_dt, relaxationupdate, val_vector);
          }
          if (apply_gamma_relaxation) {
            // bisection
            double tol = 1e-12;
            double left, right, val_left, val_right, val;
            for (double offset = 0.01; offset < 0.22; offset += 0.1) {
              left = gamma - offset;
              right = gamma + offset;
              val_left = r_gamma(left,
                                 local_alpha,
                                 local_d,
                                 grid_view,
                                 merged_quads,
                                 basis_vals,
                                 actual_dt,
                                 relaxationupdate,
                                 val_vector);
              val_right = r_gamma(right,
                                  local_alpha,
                                  local_d,
                                  grid_view,
                                  merged_quads,
                                  basis_vals,
                                  actual_dt,
                                  relaxationupdate,
                                  val_vector);
              const double val_gamma = r_gamma(gamma,
                                               local_alpha,
                                               local_d,
                                               grid_view,
                                               merged_quads,
                                               basis_vals,
                                               actual_dt,
                                               relaxationupdate,
                                               val_vector);
              if (XT::Common::FloatCmp::eq(val_gamma, 0., tol, tol)) {
                val = val_gamma;
                break;
              }
              if (XT::Common::FloatCmp::eq(val_left, 0., tol, tol)) {
                val = val_left;
                gamma = left;
                break;
              }
              if (XT::Common::FloatCmp::eq(val_right, 0., tol, tol)) {
                val = val_right;
                gamma = right;
                break;
              }
              if (val_left * val_gamma < 0.) {
                right = gamma;
                val_right = val_gamma;
                break;
              }
              if (val_right * val_gamma < 0.) {
                left = gamma;
                val_left = val_gamma;
                break;
              }
            }
            if (!XT::Common::FloatCmp::eq(val, 0., tol, tol) && val_left * val_right > 0.
                && XT::Common::FloatCmp::gt(std::max(std::abs(val_left), std::abs(val_right)), 0., tol, tol)) {
              mixed_error = 1e10;
              time_step_scale_factor = 0.5;
              write_gamma(local_alpha,
                          local_d,
                          grid_view,
                          merged_quads,
                          basis_vals,
                          actual_dt,
                          relaxationupdate,
                          val_vector,
                          "error1_gamma_");
              std::cout << "Timestep reduced due non-convex/concave r_gamma!" << std::endl;
              continue;
            }
            double pos_gamma = val_left > 0 ? left : right;
            double neg_gamma = val_left > 0 ? right : left;
            size_t num_its = 0;
            while (XT::Common::FloatCmp::ne(val, 0., tol, tol) && num_its < 300) {
              gamma = (pos_gamma + neg_gamma) / 2.;
              val = r_gamma(gamma,
                            local_alpha,
                            local_d,
                            grid_view,
                            merged_quads,
                            basis_vals,
                            actual_dt,
                            relaxationupdate,
                            val_vector);
              if (val < 0)
                neg_gamma = gamma;
              else
                pos_gamma = gamma;
              ++num_its;
              // if (num_its == 100)
              //   tol *= 10;
              // if (num_its == 200)
              //   tol *= 10;
            }
            if (num_its == 300) {
              mixed_error = 1e10;
              time_step_scale_factor = 0.5;
              std::cout << "Timestep reduced due non-convex/concave r_gamma!" << std::endl;
              write_gamma(local_alpha,
                          local_d,
                          grid_view,
                          merged_quads,
                          basis_vals,
                          actual_dt,
                          relaxationupdate,
                          val_vector,
                          "error2_gamma_");
              continue;
            }
            if (num_its != 0) {
              const std::string gamma_filename = prefix + "_gammas.txt";
              std::ofstream gamma_file(gamma_filename, std::ios_base::app);
              gamma_file << XT::Common::to_string(gamma, 15) << " " << num_its << " " << XT::Common::to_string(val, 15)
                         << std::endl;
              gamma_file.close();
            }
          } // apply_gamma_relaxation
        } // if (mixed_error < 1)
      } // if (!skip_error_computation)
    } // while (mixed_error > 1.)

    // compute alpha_np1_gamma
    if (apply_gamma_relaxation) {
      alpha_np1_->dofs().vector() = alpha_n.dofs().vector();
      actual_dt *= gamma;
      for (size_t ii = 0; ii < num_stages_ - 1; ++ii)
        alpha_np1_->dofs().vector().axpy(actual_dt * r_ * b_1_[ii], stages_k_[ii]->dofs().vector());
    }

    alpha_n.dofs().vector() = alpha_np1_->dofs().vector();
    this->dts_.push_back(actual_dt);

    if (!last_stage_of_previous_step_)
      last_stage_of_previous_step_ = alpha_n.copy_as_discrete_function();
    last_stage_of_previous_step_->dofs().vector() = stages_k_[num_stages_ - 1]->dofs().vector();
    last_relaxation_of_previous_step_ = relaxationupdate_stages_[num_stages_ - 1];
    last_gamma_ = gamma;

    t += actual_dt;

    write_entropy(local_alpha, grid_view, relaxationupdate, merged_quads, basis_vals, t, actual_dt, val_vector, prefix);

    return actual_dt * time_step_scale_factor;
  } // ... step(...)

private:
  bool check_for_nan(const RangeFieldType* vec1, const RangeFieldType* vec2, const size_t num_dofs)
  {
    for (size_t ii = 0; ii < num_dofs; ++ii)
      if (std::isnan(vec1[ii] + vec2[ii]))
        return true;
    return false;
  }

  template <class MergedQuadratureType>
  std::vector<XT::Common::FieldVector<double, dimRange>> get_basis_vals(const MergedQuadratureType& merged_quads)
  {
    const auto& basis_functions = entropy_flux_.basis_functions();
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
  double r_gamma(const double gamma,
                 LocalDiscreteFunctionType& local_alpha,
                 LocalDiscreteFunctionType& local_d,
                 const GV& grid_view,
                 const MergedQuadratureType& merged_quads,
                 const BasisValsType& basis_vals,
                 const double actual_dt,
                 const double relaxationupdate,
                 std::vector<double>& vals)
  {
    vals.clear();
    XT::Common::FieldVector<RangeFieldType, dimRange> local_alpha_vec, local_alpha_gamma_vec;
    for (auto&& entity : Dune::elements(grid_view)) {
      local_alpha->bind(entity);
      local_d->bind(entity);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        local_alpha_vec[jj] = local_alpha->dofs().get_entry(jj);
        local_alpha_gamma_vec[jj] = local_alpha_vec[jj] + gamma * local_d->dofs().get_entry(jj);
      }
      size_t kk = 0;
      for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it, ++kk) {
        const auto& quad_point = *it;
        const auto alpha_n_b = local_alpha_vec * basis_vals[kk];
        const auto alpha_gamma_b = local_alpha_gamma_vec * basis_vals[kk];
        vals.push_back(std::exp(alpha_gamma_b) * (alpha_gamma_b - 1) * quad_point.weight());
        vals.push_back(-std::exp(alpha_n_b) * (alpha_n_b - 1) * quad_point.weight());
      } // quad_points
    } // entities
    vals.push_back(-gamma * actual_dt * relaxationupdate);
    return compensated_sum(vals);
  }

  template <class GV, class LocalDiscreteFunctionType, class MergedQuadratureType, class BasisValsType>
  double compute_entropy(LocalDiscreteFunctionType& local_alpha,
                         const GV& grid_view,
                         const MergedQuadratureType& merged_quads,
                         const BasisValsType& basis_vals,
                         std::vector<double>& vals)
  {
    vals.clear();
    XT::Common::FieldVector<RangeFieldType, dimRange> local_alpha_vec;
    for (auto&& entity : Dune::elements(grid_view)) {
      local_alpha->bind(entity);
      for (size_t jj = 0; jj < dimRange; ++jj)
        local_alpha_vec[jj] = local_alpha->dofs().get_entry(jj);
      size_t kk = 0;
      for (auto it = merged_quads.begin(); it != merged_quads.end(); ++it, ++kk) {
        const auto& quad_point = *it;
        const auto alpha_n_b = local_alpha_vec * basis_vals[kk];
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

  template <class GV, class LocalDiscreteFunctionType, class MergedQuadratureType, class BasisValsType>
  void write_gamma(LocalDiscreteFunctionType& local_alpha,
                   LocalDiscreteFunctionType& local_d,
                   const GV& grid_view,
                   const MergedQuadratureType& merged_quads,
                   const BasisValsType& basis_vals,
                   const double actual_dt,
                   const double relaxationupdate,
                   std::vector<double>& val_vector,
                   std::string prefix = "gamma_")
  {
    static size_t ll = 0;
    const std::string gamma_filename = prefix + XT::Common::to_string(ll) + ".txt";
    std::ofstream gamma_file(gamma_filename);
    for (double gamma2 = 0.; gamma2 < 1.209; gamma2 += 0.01) {
      double val = r_gamma(
          gamma2, local_alpha, local_d, grid_view, merged_quads, basis_vals, actual_dt, relaxationupdate, val_vector);
      gamma_file << XT::Common::to_string(gamma2) << " " << XT::Common::to_string(val, 15) << std::endl;
    } // gamma
    gamma_file.close();
    ++ll;
  }

  template <class GV, class LocalDiscreteFunctionType, class MergedQuadratureType, class BasisValsType>
  void write_entropy(LocalDiscreteFunctionType& local_alpha,
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
    const double entropy = compute_entropy(local_alpha, grid_view, merged_quads, basis_vals, val_vector);
    const double diff =
        XT::Common::is_zero(t) ? 0. : compensated_sum({entropy, -1. * last_entropy_, -dt * relaxationupdate});
    entropy_file << XT::Common::to_string(t) << " " << XT::Common::to_string(entropy, 15) << " "
                 << XT::Common::to_string(diff, 15) << " " << XT::Common::to_string(std::abs(diff), 15) << std::endl;
    entropy_file.close();
    last_entropy_ = entropy;
  }

  const MinDensitySetterType min_density_setter_;
  const OperatorType& op_;
  const EntropyFluxType& entropy_flux_;
  const RangeFieldType r_;
  const RangeFieldType atol_;
  const RangeFieldType rtol_;
  const RangeFieldType scale_factor_min_;
  const RangeFieldType scale_factor_max_;
  std::unique_ptr<DiscreteFunctionType> alpha_tmp_;
  std::unique_ptr<DiscreteFunctionType> alpha_np1_;
  const MatrixType A_;
  const VectorType b_1_;
  const VectorType b_2_;
  const VectorType c_;
  const VectorType b_diff_;
  std::vector<std::unique_ptr<DiscreteFunctionType>> stages_k_;
  const size_t num_stages_;
  std::unique_ptr<DiscreteFunctionType> last_stage_of_previous_step_;
  double last_relaxation_of_previous_step_;
  bool first_same_as_last_;
  std::vector<double> relaxationupdate_stages_;
  double last_gamma_;
  double last_entropy_;
  XT::Common::Parameter op_param_;
}; // class KineticAdaptiveRungeKuttaTimeStepper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_KINETIC_ADAPTIVE_RUNGEKUTTA_HH
