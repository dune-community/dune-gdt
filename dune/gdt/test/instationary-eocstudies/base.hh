// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_BASE_HH
#define DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_BASE_HH

#include <cmath>
#include <functional>
#include <memory>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/test/gtest/gtest.h>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/container/vector-array/list.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/base/sliced.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/functionals/localizable-functional.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/abs.hh>
#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/constant.hh>
#include <dune/gdt/operators/identity.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/prolongations.hh>
#include <dune/gdt/spaces/bochner.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class GridView, size_t m_ = 1, XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class InstationaryEocStudy
  : public XT::Common::ConvergenceStudy
  , public ::testing::Test
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

protected:
  using GV = GridView;
  using G = typename GridView::Grid;
  static const constexpr size_t m = m_;
  using D = double;
  static const constexpr size_t d = G::dimension;
  using R = double;
  using E = XT::Grid::extract_entity_t<GV>;
  using DomainType = XT::Common::FieldVector<D, d>;
  using RangeType = XT::Common::FieldVector<D, m>;
  using GP = XT::Grid::GridProvider<G>;
  using S = SpaceInterface<GV, m>;
  using M = typename XT::LA::Container<R, la>::MatrixType;
  using V = XT::LA::vector_t<M>;
  using DF = DiscreteFunction<V, GV, m>;
  using BS = BochnerSpace<GV, m>;
  using O = OperatorInterface<M, GV, m>;

public:
  InstationaryEocStudy(
      const double T_end,
      const std::string timestepping,
      std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualizer =
          [](const auto& /*solution*/, const auto& /*prefix*/) { /*no visualization by default*/ },
      const size_t num_refinements = DXTC_CONFIG_GET("num_refinements", 3),
      const size_t num_additional_refinements_for_reference =
          DXTC_CONFIG_GET("num_additional_refinements_for_reference", 2))
    : T_end_(T_end)
    , timestepping_(timestepping)
    , num_refinements_(num_refinements)
    , num_additional_refinements_for_reference_(num_additional_refinements_for_reference)
    , visualize_(visualizer)
    , current_refinement_(std::numeric_limits<size_t>::max())
    , current_data_()
    , current_grid_(nullptr)
    , current_space_(nullptr)
    , reference_grid_(nullptr)
    , reference_space_(nullptr)
    , current_solution_on_reference_grid_(nullptr)
    , reference_solution_on_reference_grid_(nullptr)
  {}

  size_t num_refinements() const override
  {
    return num_refinements_;
  }

  std::vector<std::string> targets() const override
  {
    return {"h"};
  }

  std::vector<std::string> norms() const override
  {
    // We currently support the following temporal norms:
    //   L_infty
    //   L_2
    // and the following spatial norms:
    //   L_1
    //   L_2
    return {"L_infty/L_2"};
  }

  std::vector<std::pair<std::string, std::string>> estimates() const override
  {
    return {};
  }

  std::vector<std::string> quantities() const override
  {
    std::vector<std::string> ret = {"time to solution (s)", "rel mass conserv error", "num timesteps"};
    if (this->adaptive_timestepping()) {
      ret.push_back("min dt");
      ret.push_back("max dt");
    } else {
      ret.push_back("CFL");
    }
    return ret;
  } // ... quantities(...)

  std::string discretization_info_title() const override
  {
    return " |grid| |   #DoFs";
  }

protected:
  virtual GP make_initial_grid() = 0;

  virtual std::unique_ptr<S> make_space(const GP& current_grid) = 0;

public:
  std::string discretization_info(const size_t refinement_level) override
  {
    if (current_refinement_ != refinement_level) {
      // clear the current state
      current_grid_.reset();
      current_space_.reset();
      current_solution_on_reference_grid_.reset();
      // compute on this refinement
      current_grid_ = std::make_unique<GP>(make_initial_grid());
      for (size_t ref = 0; ref < refinement_level; ++ref)
        current_grid_->global_refine(DGFGridInfo<G>::refineStepsForHalf());
      current_space_ = make_space(*current_grid_);
      double grid_size = 0;
      double grid_width = 0.;
      for (auto&& grid_element : elements(current_space_->grid_view())) {
        grid_size += 1;
        grid_width = std::max(grid_width, XT::Grid::entity_diameter(grid_element));
      }
      current_refinement_ = refinement_level;
      current_data_.clear();
      current_data_["target"]["h"] = grid_width;
      current_data_["quantity"]["num_grid_elements"] = grid_size;
      current_data_["quantity"]["num_dofs"] = current_space_->mapper().size();
    }
    const auto lfill_nicely = [&](const auto& number, const auto& len) {
      std::string ret;
      using XT::Common::to_string;
      if (std::log10(number) < len)
        ret = this->lfill(to_string(number), len);
      else {
        std::stringstream ss;
        ss << std::setprecision(2) << std::scientific << number;
        ret = this->lfill(ss.str(), len);
      }
      return ret;
    };
    return lfill_nicely(current_data_["quantity"]["num_grid_elements"], 7) + " | "
           + lfill_nicely(current_data_["quantity"]["num_dofs"], 7);
  } // ... discretization_info(...)

protected:
  virtual std::unique_ptr<O> make_lhs_operator(const S& space) = 0;

  virtual XT::LA::ListVectorArray<V> solve(const S& space, const double T_end) = 0;

  virtual void compute_reference_solution()
  {
    if (reference_solution_on_reference_grid_)
      return;
    reference_grid_ = std::make_unique<GP>(make_initial_grid());
    for (size_t ref = 0; ref < num_refinements_ + num_additional_refinements_for_reference_; ++ref)
      reference_grid_->global_refine(DGFGridInfo<G>::refineStepsForHalf());
    reference_space_ = make_space(*reference_grid_);
    reference_solution_on_reference_grid_ =
        std::make_unique<XT::LA::ListVectorArray<V>>(solve(*reference_space_, T_end_));
    // visualize
    const BS reference_bochner_space(*reference_space_,
                                     time_points_from_vector_array(*reference_solution_on_reference_grid_));
    visualize_(make_discrete_bochner_function(reference_bochner_space, *reference_solution_on_reference_grid_),
               "reference_solution_on_refinement_"
                   + XT::Common::to_string(num_refinements_ + num_additional_refinements_for_reference_));
  } // ... compute_reference_solution(...)

public:
  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level,
          const std::vector<std::string>& actual_norms,
          const std::vector<std::pair<std::string, std::string>>& actual_estimates,
          const std::vector<std::string>& actual_quantities) override
  {
    auto& self = *this;
    if (current_refinement_ != refinement_level)
      self.discretization_info(refinement_level);
    DUNE_THROW_IF(!current_space_, InvalidStateException, "");
    // compute current solution
    const auto& current_space = *current_space_;
    Timer timer;
    auto solution_on_current_grid = solve(current_space, T_end_);
    const double time_to_solution = timer.elapsed();
    // only set time if this did not happen in solve()
    if (current_data_["quantity"].count("time to solution (s)") == 0)
      current_data_["quantity"]["time to solution (s)"] = time_to_solution;
    // visualize
    const BS current_bochner_space(current_space, time_points_from_vector_array(solution_on_current_grid));
    visualize_(make_discrete_bochner_function(current_bochner_space, solution_on_current_grid),
               "solution_on_refinement_" + XT::Common::to_string(refinement_level));
    const auto coarse_solution = make_discrete_bochner_function(current_bochner_space, solution_on_current_grid);
    // determine wether we need a reference solution
    // compute statistics
    auto norms_to_compute = actual_norms;
    // - norms
    if (!norms_to_compute.empty()) {
      auto current_data_backup = current_data_;
      self.compute_reference_solution();
      current_data_ = current_data_backup;
      DUNE_THROW_IF(!reference_space_, InvalidStateException, "");
      const auto& reference_space = *reference_space_;
      DUNE_THROW_IF(!reference_solution_on_reference_grid_, InvalidStateException, "");
      auto& u = *reference_solution_on_reference_grid_;
      // prolong
      const BS reference_bochner_space(reference_space,
                                       time_points_from_vector_array(*reference_solution_on_reference_grid_));
      current_solution_on_reference_grid_ = std::make_unique<XT::LA::ListVectorArray<V>>(
          prolong<V>(coarse_solution, reference_bochner_space).dof_vectors());
      auto& u_h = *current_solution_on_reference_grid_;
      while (!norms_to_compute.empty()) {
        const auto norm_id = norms_to_compute.back();
        const auto components = XT::Common::tokenize(norm_id, "/");
        norms_to_compute.pop_back();
        // compute Bochner norm
        DUNE_THROW_IF(components.size() != 2,
                      XT::Common::Exceptions::wrong_input_given,
                      "I do not know how to compute the following norm: " << norm_id);
        const auto temporal_norm_id = components[0];
        const auto spatial_norm_id = components[1];
        // - spatial component
        std::function<double(const DF&)> spatial_norm;
        if (spatial_norm_id == "L_1") {
          spatial_norm = [&](const DF& func) {
            auto localizable_functional = make_localizable_functional(reference_space.grid_view(), func);
            localizable_functional.append(LocalElementIntegralFunctional<E, m>(LocalElementAbsIntegrand<E, m>()));
            localizable_functional.assemble();
            return localizable_functional.result();
          };
        } else if (spatial_norm_id == "L_2") {
          spatial_norm = [&](const DF& func) {
            auto localizable_product = make_localizable_bilinear_form(reference_space.grid_view(), func, func);
            localizable_product.append(LocalElementIntegralBilinearForm<E, m>(LocalElementProductIntegrand<E, m>()));
            localizable_product.assemble();
            return std::sqrt(localizable_product.result());
          };
        } else
          DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                     "I do not know how to compute the spatial norm '" << spatial_norm_id << "'!");
        // - temporal component
        if (temporal_norm_id == "L_infty") {
          double result = 0;
          for (size_t ii = 0; ii < u.length(); ++ii)
            result = std::max(result,
                              spatial_norm(make_discrete_function(reference_space, u[ii].vector() - u_h[ii].vector())));
          current_data_["norm"][norm_id] = result;
        } else if (temporal_norm_id == "L_2") {
          const XT::Functions::GenericFunction<1> spatial_norm_function(
              1, [&](const auto& time, const auto& /*param*/) {
                const auto u_t = make_discrete_bochner_function(reference_bochner_space, u).evaluate(time);
                const auto u_h_t = make_discrete_bochner_function(reference_bochner_space, u_h).evaluate(time);
                return spatial_norm(
                    make_discrete_function(reference_space, u_t.dofs().vector() - u_h_t.dofs().vector()));
              });
          auto temporal_grid_view = reference_bochner_space.temporal_space().grid_view();
          using TE = XT::Grid::extract_entity_t<decltype(temporal_grid_view)>;
          const auto& spatial_norm_as_temporal_grid_function = spatial_norm_function.template as_grid_function<TE>();
          auto localizable_temporal_product = make_localizable_bilinear_form(
              temporal_grid_view, spatial_norm_as_temporal_grid_function, spatial_norm_as_temporal_grid_function);
          localizable_temporal_product.append(LocalElementIntegralBilinearForm<TE>(LocalElementProductIntegrand<TE>()));
          localizable_temporal_product.assemble();
          current_data_["norm"][norm_id] = localizable_temporal_product.result();
        } else
          DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                     "I do not know how to compute the temporal norm '" << temporal_norm_id << "'!");
      }
    }
    DUNE_THROW_IF(!norms_to_compute.empty(),
                  XT::Common::Exceptions::wrong_input_given,
                  "I did not know how to compute the following norms: " << norms_to_compute);
    // - estimates
    auto estimats_to_compute = actual_estimates;
    DUNE_THROW_IF(!estimats_to_compute.empty(),
                  XT::Common::Exceptions::wrong_input_given,
                  "I did not know how to compute the following estimates: " << estimats_to_compute);
    // - quantities
    auto quantities_to_compute = actual_quantities;
    while (!quantities_to_compute.empty()) {
      const auto id = quantities_to_compute.back();
      quantities_to_compute.pop_back();
      if (id == "time to solution (s)") {
        DUNE_THROW_IF(current_data_["quantity"].find(id) == current_data_["quantity"].end(),
                      InvalidStateException,
                      "Could not find id " << id << "in current_data_ map");
      } else if (id == "rel mass conserv error") {
        const auto compute_masses = [&](const auto& vec) {
          auto func = make_discrete_function(current_space, vec);
          FieldVector<double, m> results(0);
          for (size_t ss = 0; ss < m; ++ss) {
            auto component_function = XT::Functions::make_sliced_function<1>(func, {ss});
            auto l1_functional = make_localizable_functional(current_space.grid_view(), component_function);
            l1_functional.append(LocalElementIntegralFunctional<E>(LocalElementIdentityIntegrand<E>()));
            l1_functional.assemble(/*use_tbb=*/true);
            results[ss] = l1_functional.result();
          }
          return results;
        }; // ... compute_masses(...)
        const auto initial_masses = compute_masses(solution_on_current_grid[0].vector());
        auto relative_mass_conservation_errors = initial_masses;
        relative_mass_conservation_errors *= 0;
        for (size_t ii = 1; ii < solution_on_current_grid.length(); ++ii) {
          const auto current_masses = compute_masses(solution_on_current_grid[ii].vector());
          for (size_t ss = 0; ss < m; ++ss)
            relative_mass_conservation_errors[ss] = std::max(relative_mass_conservation_errors[ss],
                                                             std::abs(initial_masses[ss] - current_masses[ss])
                                                                 / (initial_masses[ss] > 0. ? initial_masses[ss] : 1.));
        }
        current_data_["quantity"][id] = relative_mass_conservation_errors.infinity_norm();
      } else if (id == "num timesteps") {
        current_data_["quantity"][id] = solution_on_current_grid.vectors().size();
      } else if (id == "CFL") {
        DUNE_THROW_IF(this->adaptive_timestepping(), InvalidStateException, "");
        const auto time_points = this->time_points_from_vector_array(solution_on_current_grid);
        const auto dt = time_points.at(1) - time_points.at(0);
        current_data_["quantity"][id] = dt / this->extract(current_data_, "quantity", "explicit_fv_dt");
      } else if (id == "dt") {
        DUNE_THROW_IF(this->adaptive_timestepping(), InvalidStateException, "");
        const auto time_points = this->time_points_from_vector_array(solution_on_current_grid);
        current_data_["quantity"][id] = time_points.at(1) - time_points.at(0);
      } else if (id == "min dt") {
        double min_dt = std::numeric_limits<double>::max();
        const auto time_points = this->time_points_from_vector_array(solution_on_current_grid);
        for (size_t ii = 1; ii < time_points.size(); ++ii) {
          auto dt = time_points[ii] - time_points[ii - 1];
          min_dt = std::min(min_dt, dt);
        }
        current_data_["quantity"][id] = min_dt;
      } else if (id == "max dt") {
        double max_dt = std::numeric_limits<double>::min();
        const auto time_points = this->time_points_from_vector_array(solution_on_current_grid);
        for (size_t ii = 1; ii < time_points.size(); ++ii) {
          auto dt = time_points[ii] - time_points[ii - 1];
          max_dt = std::max(max_dt, dt);
        }
        current_data_["quantity"][id] = max_dt;
      } else
        DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                   "I do not know how to compute the quantity '" << id << "'!");
    }
    DUNE_THROW_IF(!quantities_to_compute.empty(),
                  XT::Common::Exceptions::wrong_input_given,
                  "I did not know how to compute the following quantities: " << quantities_to_compute);
    return current_data_;
  } // ... compute_on_current_refinement(...)

protected:
  static std::vector<double> time_points_from_vector_array(const XT::LA::ListVectorArray<V>& va)
  {
    std::vector<double> time_points;
    time_points.reserve(va.length());
    for (const auto& note : va.notes())
      time_points.emplace_back(note.get("_t").at(0));
    return time_points;
  }

  bool implicit_timestepping() const
  {
    auto ts = XT::Common::tokenize(timestepping_, "/");
    DUNE_THROW_IF(ts.size() != 2, XT::Common::Exceptions::wrong_input_given, ts);
    DUNE_THROW_IF(ts[0] != "explicit" && ts[0] != "implicit", XT::Common::Exceptions::wrong_input_given, ts[0]);
    return ts[0] == "implicit";
  }

  bool adaptive_timestepping() const
  {
    auto ts = XT::Common::tokenize(timestepping_, "/");
    DUNE_THROW_IF(ts.size() != 2, XT::Common::Exceptions::wrong_input_given, ts);
    DUNE_THROW_IF(ts[1] != "fixed" && ts[1] != "adaptive", XT::Common::Exceptions::wrong_input_given, ts[1]);
    return ts[1] == "adaptive";
  }

  double T_end_;
  const std::string timestepping_;
  size_t num_refinements_;
  size_t num_additional_refinements_for_reference_;
  const std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualize_;
  size_t current_refinement_;
  std::map<std::string, std::map<std::string, double>> current_data_;
  std::unique_ptr<GP> current_grid_;
  std::unique_ptr<S> current_space_;
  std::unique_ptr<GP> reference_grid_;
  std::unique_ptr<S> reference_space_;
  std::unique_ptr<XT::LA::ListVectorArray<V>> current_solution_on_reference_grid_;
  std::unique_ptr<XT::LA::ListVectorArray<V>> reference_solution_on_reference_grid_;
}; // struct InstationaryEocStudy


template <class V, class GV, size_t m, class M>
XT::LA::ListVectorArray<V> solve_instationary_system_explicit_euler(const DiscreteFunction<V, GV, m>& initial_values,
                                                                    const OperatorInterface<M, GV, m>& spatial_op,
                                                                    const double T_end,
                                                                    const double dt)
{
  // initial values
  XT::LA::ListVectorArray<V> solution(
      spatial_op.source_space().mapper().size(), /*length=*/0, /*reserve=*/std::ceil(T_end / (dt)));
  solution.append(initial_values.dofs().vector(), {"_t", 0.});
  // timestepping
  double time = 0.;
  while (time < T_end + dt) {
    const auto& u_n = solution.back().vector();
    auto u_n_plus_one = u_n - spatial_op.apply(u_n, {{"_t", {time}}, {"_dt", {dt}}}) * dt;
    time += dt;
    solution.append(std::move(u_n_plus_one), {"_t", time});
  }
  return solution;
} // ... solve_instationary_system_explicit_euler(...)


template <class V, class GV, size_t m, class M>
XT::LA::ListVectorArray<V> solve_instationary_system_implicit_euler(const DiscreteFunction<V, GV, m>& initial_values,
                                                                    const OperatorInterface<M, GV, m>& spatial_op,
                                                                    const double T_end,
                                                                    const double dt)
{
  // some preparations
  auto id = make_identity_operator(spatial_op);
  V zero(spatial_op.range_space().mapper().size(), 0.);
  auto logger = XT::Common::TimedLogger().get("gdt.test.solve_instationary_system_implicit_euler");
  // initial values
  XT::LA::ListVectorArray<V> solution(
      spatial_op.source_space().mapper().size(), /*length=*/0, /*reserve=*/std::ceil(T_end / (dt)));
  solution.append(initial_values.dofs().vector(), {"_t", 0.});
  // timestepping
  double time = 0.;
  while (time < T_end + dt) {
    logger.debug() << "time = " << time << ": stepping with dt=" << dt << "..." << std::endl;
    time += dt;
    const auto& u_n = solution.back().vector();
    auto residual_op = (id - u_n) / dt + spatial_op;
    auto u_n_plus_one = u_n.copy();
    residual_op.apply_inverse(zero,
                              u_n_plus_one,
                              DXTC_CONFIG.has_sub("solve_instationary_system_implicit_euler.apply_inverse")
                                  ? DXTC_CONFIG.sub("solve_instationary_system_implicit_euler.apply_inverse")
                                  : residual_op.invert_options(residual_op.invert_options().at(0)));
    solution.append(std::move(u_n_plus_one), {"_t", time});
  }
  return solution;
} // ... solve_instationary_system_implicit_euler(...)


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_BASE_HH
