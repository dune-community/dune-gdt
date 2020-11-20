// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_BASE_HH
#define DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_BASE_HH

#include <cmath>
#include <functional>
#include <memory>

#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dune/xt/common/convergence-study.hh>
#include <dune/xt/test/common.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/string.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/test/common.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/functionals/localizable-functional.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/local/integrands/abs.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/norms.hh>
#include <dune/gdt/operators/constant.hh>
#include <dune/gdt/operators/identity.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/bilinear-form.hh>
#include <dune/gdt/operators/matrix.hh>
#include <dune/gdt/prolongations.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class GridView, size_t m_ = 1, XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class StationaryEocStudy
  : public XT::Common::ConvergenceStudy
  , public ::testing::Test
{
  static_assert(XT::Grid::is_view<GridView>::value, "");

protected:
  using GV = GridView;
  using G = typename GridView::Grid;
  static constexpr size_t m = m_;
  using D = double;
  static constexpr size_t d = G::dimension;
  using R = double;
  using I = XT::Grid::extract_intersection_t<GV>;
  using DomainType = XT::Common::FieldVector<D, d>;
  using RangeType = XT::Common::FieldVector<D, m>;
  using GP = XT::Grid::GridProvider<G>;
  using S = SpaceInterface<GV, m>;
  using M = typename XT::LA::Container<R, la>::MatrixType;
  using V = XT::LA::vector_t<M>;
  using DF = DiscreteFunction<V, GV, m>;
  using O = OperatorInterface<M, GV, m>;

public:
  using E = XT::Grid::extract_entity_t<GV>;

  StationaryEocStudy(
      std::function<void(const DF&, const std::string&)> visualizer =
          [](const auto& solution, const auto& prefix) {
            if (DXTC_TEST_CONFIG_GET("setup.visualize", false))
              solution.visualize(prefix);
          },
      const size_t num_refinements = DXTC_TEST_CONFIG_GET("setup.num_refinements", 3),
      const size_t num_additional_refinements_for_reference =
          DXTC_TEST_CONFIG_GET("setup.num_additional_refinements_for_reference", 2))
    : num_refinements_(num_refinements)
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
    // We currently support the following norms: L_1, L_2, H_1_semi
    return {"L_2", "H_1_semi"};
  }

  std::vector<std::pair<std::string, std::string>> estimates() const override
  {
    return {};
  }

  std::vector<std::string> quantities() const override
  {
    return {"time to solution (s)"};
  }

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
      current_data_["quantity"]["num_dofs"] = static_cast<double>(current_space_->mapper().size());
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
  virtual std::unique_ptr<O> make_residual_operator(const S& space) = 0;

  virtual V solve(const S& space)
  {
    auto op = make_residual_operator(space);
    V zero(op->range_space().mapper().size(), 0.);
    return op->apply_inverse(zero);
  }

  virtual void compute_reference_solution()
  {
    if (reference_solution_on_reference_grid_)
      return;
    reference_grid_ = std::make_unique<GP>(make_initial_grid());
    for (size_t ref = 0; ref < num_refinements_ + num_additional_refinements_for_reference_; ++ref)
      reference_grid_->global_refine(DGFGridInfo<G>::refineStepsForHalf());
    reference_space_ = make_space(*reference_grid_);
    reference_solution_on_reference_grid_ = std::make_unique<V>(std::move(solve(*reference_space_)));
    // visualize
    visualize_(make_discrete_function(*reference_space_, *reference_solution_on_reference_grid_),
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
    auto solution_on_current_grid = solve(current_space);
    // only set time if this did not happen in solve()
    if (current_data_["quantity"].count("time to solution (s)") == 0)
      current_data_["quantity"]["time to solution (s)"] = timer.elapsed();
    // visualize
    visualize_(make_discrete_function(current_space, solution_on_current_grid),
               "solution_on_refinement_" + XT::Common::to_string(refinement_level));
    const auto coarse_solution = make_discrete_function(current_space, solution_on_current_grid);
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
      current_solution_on_reference_grid_ =
          std::make_unique<V>(prolong<V>(coarse_solution, reference_space).dofs().vector());
      // compute norm
      auto& u_h = *current_solution_on_reference_grid_;
      while (!norms_to_compute.empty()) {
        const auto spatial_norm_id = norms_to_compute.back();
        norms_to_compute.pop_back();
        std::function<double(const DF&)> spatial_norm;
        if (spatial_norm_id == "L_1") {
          spatial_norm = [&](const DF& func) {
            auto localizable_functional = make_localizable_functional(reference_space.grid_view(), func);
            localizable_functional.append(LocalElementIntegralFunctional<E, m>(
                LocalElementAbsIntegrand<E, m>(), DXTC_TEST_CONFIG_GET("setup.over_integrate", 3)));
            localizable_functional.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
            return localizable_functional.result();
          };
        } else if (spatial_norm_id == "L_2") {
          spatial_norm = [&](const DF& func) {
            return l2_norm(reference_space.grid_view(), func, DXTC_TEST_CONFIG_GET("setup.over_integrate", 3));
          };
        } else if (spatial_norm_id == "H_1_semi") {
          spatial_norm = [&](const DF& func) {
            auto product = make_bilinear_form(reference_space.grid_view(), func, func);
            product += LocalElementIntegralBilinearForm<E, m>(LocalLaplaceIntegrand<E, m>(),
                                                              DXTC_TEST_CONFIG_GET("setup.over_integrate", 3));
            return std::sqrt(product.apply2());
          };
        } else
          DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                     "I do not know how to compute the norm '" << spatial_norm_id << "'!");
        current_data_["norm"][spatial_norm_id] = spatial_norm(make_discrete_function(reference_space, u - u_h));
      }
      DUNE_THROW_IF(!norms_to_compute.empty(),
                    XT::Common::Exceptions::wrong_input_given,
                    "I did not know how to compute the following norms: " << norms_to_compute);
    }
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
  size_t num_refinements_;
  size_t num_additional_refinements_for_reference_;
  const std::function<void(const DF&, const std::string&)> visualize_;
  size_t current_refinement_;
  std::map<std::string, std::map<std::string, double>> current_data_;
  std::unique_ptr<GP> current_grid_;
  std::unique_ptr<S> current_space_;
  std::unique_ptr<GP> reference_grid_;
  std::unique_ptr<S> reference_space_;
  std::unique_ptr<V> current_solution_on_reference_grid_;
  std::unique_ptr<V> reference_solution_on_reference_grid_;
}; // struct StationaryEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_EOCSTUDIES_BASE_HH
