// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/instationary-eocstudies/hyperbolic-fv.hh>

#include "base.hh"

using namespace Dune;
using namespace Dune::GDT::Test;


template <class G>
class LinearTransportExplicitFvTest : XT::Common::ConstStorageProvider<LinearTransportProblem<G>>,
                                      public InstationaryHyperbolicFiniteVolumeEocStudy<G, 1>
{
  using Problem = XT::Common::ConstStorageProvider<LinearTransportProblem<G>>;
  using BaseType = InstationaryHyperbolicFiniteVolumeEocStudy<G, 1>;

protected:
  using BaseType::d;
  using typename BaseType::F;
  using typename BaseType::DF;
  using typename BaseType::GP;
  using typename BaseType::S;
  using typename BaseType::BS;
  using typename BaseType::V;

public:
  LinearTransportExplicitFvTest(const size_t num_refinements = 3,
                                const size_t num_additional_refinements_for_reference = 1) // Not exact, but enough.
      : Problem(new LinearTransportProblem<G>()),
        BaseType(Problem::access().T_end,
                 num_refinements,
                 num_additional_refinements_for_reference,
                 [&](const auto& solution, const auto& prefix) {
                   for (size_t ii = 0; ii < this->visualization_steps_; ++ii) {
                     const double time = ii * (this->T_end_ / this->visualization_steps_);
                     const auto u_t = solution.evaluate(time);
                     make_discrete_function(solution.space().spatial_space(), u_t)
                         .visualize(prefix + "_solution_" + XT::Common::to_string(ii));
                   }
                 }),
        visualization_steps_(0)
  {
  }

  std::vector<std::string> targets() const override final
  {
    if (d == 1) // in 1d dt depends linearly on h, so no need to pollute the EOC table with dt-related values
      return {"h"};
    else
      return BaseType::targets();
  }

  double estimate_dt(const S& space) override final
  {
    if (d == 1) { // in 1d we know that dt = h is a good choice
      double grid_width = 0.;
      for (auto&& grid_element : elements(space.grid_view())) {
        grid_width = std::max(grid_width, XT::Grid::entity_diameter(grid_element));
      }
      return grid_width;
    } else
      return BaseType::estimate_dt(space);
  } // ... estimate_dt(...)

protected:
  const F& flux() const override final
  {
    return Problem::access().flux;
  }

  DF make_initial_values(const S& space) override final
  {
    return Problem::access().template make_exact_solution__periodic_boundaries<V>(space, 0.);
  }

  GP make_initial_grid() override final
  {
    return Problem::access().make_initial_grid();
  }

  XT::LA::ListVectorArray<V> solve(const S& space, const double T_end, const double dt) override final
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    return solve_instationary_system_explicit_euler(u_0, *op, T_end, dt);
  }

  void compute_reference_solution() override final
  {
    auto& self = *this;
    if (self.reference_solution_on_reference_grid_)
      return;
    self.reference_grid_ = std::make_unique<GP>(make_initial_grid());
    for (size_t ref = 0; ref < self.num_refinements_ + self.num_additional_refinements_for_reference_; ++ref)
      self.reference_grid_->global_refine(DGFGridInfo<G>::refineStepsForHalf());
    self.reference_space_ = self.make_space(*self.reference_grid_);
    const auto dt = estimate_dt(*self.reference_space_);
    self.reference_solution_on_reference_grid_ = std::make_unique<XT::LA::ListVectorArray<V>>(
        self.reference_space_->mapper().size(), /*length=*/0, /*reserve=*/std::ceil(self.T_end_ / (dt)));
    double time = 0.;
    while (time < self.T_end_ + dt) {
      auto u_t = Problem::access().template make_exact_solution__periodic_boundaries<V>(*self.reference_space_, time);
      self.reference_solution_on_reference_grid_->append(u_t.dofs().vector(), {"_t", time});
      time += dt;
    }
    // visualize
    const BS reference_bochner_space(*self.reference_space_,
                                     self.time_points_from_vector_array(*self.reference_solution_on_reference_grid_));
    self.visualize_(
        make_discrete_bochner_function(reference_bochner_space, *self.reference_solution_on_reference_grid_),
        "reference_solution_on_refinement_"
            + XT::Common::to_string(self.num_refinements_ + self.num_additional_refinements_for_reference_));
  } // ... compute_reference_solution(...)

protected:
  size_t visualization_steps_;
}; // class LinearTransportExplicitFvTest


using LinearTransport1dExplicitFvTest = LinearTransportExplicitFvTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_upwind_flux)
{
  //  this->visualization_steps_ = 100; // <- something like this to visualize
  this->set_numerical_flux("upwind");
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_lax_riedrichs_flux)
{
  this->set_numerical_flux("lax_friedrichs");
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_engquist_osher_flux)
{
  this->set_numerical_flux("engquist_osher");
  this->run();
}
TEST_F(LinearTransport1dExplicitFvTest, periodic_boundaries__numerical_vijayasundaram_flux)
{
  this->set_numerical_flux("vijayasundaram");
  this->run();
}
