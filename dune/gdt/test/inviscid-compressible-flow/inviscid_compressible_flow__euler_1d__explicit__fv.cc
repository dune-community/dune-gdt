// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/functions/base/sliced.hh>

#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/test/instationary-eocstudies/hyperbolic-fv.hh>

#include "base.hh"

using namespace Dune;
using namespace Dune::GDT;
using namespace Dune::GDT::Test;


template <class G>
class InviscidCompressibleFlowEulerExplicitFvTest
    : protected XT::Common::ConstStorageProvider<InviscidCompressibleFlowEulerProblem<G>>,
      public InstationaryHyperbolicFiniteVolumeEocStudy<G, G::dimension + 2>
{
  using BaseType = InstationaryHyperbolicFiniteVolumeEocStudy<G, G::dimension + 2>;

protected:
  using Problem = XT::Common::ConstStorageProvider<InviscidCompressibleFlowEulerProblem<G>>;
  using BaseType::d;
  using BaseType::m;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::DF;
  using typename BaseType::GV;
  using typename BaseType::GP;
  using typename BaseType::S;
  using typename BaseType::BS;
  using typename BaseType::V;
  using typename BaseType::O;

public:
  InviscidCompressibleFlowEulerExplicitFvTest(const size_t num_refinements = 3,
                                              const size_t num_additional_refinements_for_reference = 2)
    : Problem(new InviscidCompressibleFlowEulerProblem<G>())
    , BaseType(Problem::access().T_end,
               num_refinements,
               num_additional_refinements_for_reference,
               [&](const auto& solution, const auto& prefix) {
                 for (size_t ii = 0; ii < this->visualization_steps_; ++ii) {
                   const double time = ii * (this->T_end_ / this->visualization_steps_);
                   const auto u_t = solution.evaluate(time);
                   Problem::access().euler_tools.visualize(
                       u_t, u_t.space().grid_view(), prefix, XT::Common::to_string(ii));
                 }
               })
    , visualization_steps_(0)
    , boundary_treatment("")
  {
  }

  std::vector<std::string> targets() const override final
  {
    if (d == 1) // in 1d dt depends linearly on h, so no need to pollute the EOC table with dt-related values
      return {"h"};
    else
      return BaseType::targets();
  }

  std::vector<std::string> quantities() const override final
  {
    return {"rel. mass conserv. error"};
  }

  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level, const std::vector<std::string>& only_these) override final
  {
    auto& self = *this;
    auto data = BaseType::compute(refinement_level, only_these);
    DUNE_THROW_IF(!self.reference_space_, InvalidStateException, "");
    const auto& space = *self.reference_space_;
    DUNE_THROW_IF(!self.current_solution_on_reference_grid_, InvalidStateException, "");
    const auto& u = *self.current_solution_on_reference_grid_;
    auto actual_quantities = self.filter(self.quantities(), only_these);
    while (!actual_quantities.empty()) {
      const auto id = actual_quantities.back();
      actual_quantities.pop_back();
      if (id == "rel. mass conserv. error") {
        const auto compute_mass = [&](const auto& vec) {
          const auto func = make_discrete_function(space, vec);
          const auto density = XT::Functions::make_sliced_function<1>(func, {0}, "density");
          auto localizable_functional = make_localizable_functional(space.grid_view(), density);
          localizable_functional.append(LocalElementIntegralFunctional<E>(LocalElementIdentityIntegrand<E>()));
          localizable_functional.assemble();
          return localizable_functional.result();
        };
        const double initial_mass = compute_mass(u[0].vector());
        double relative_mass_conservation_error = 0.;
        for (size_t ii = 1; ii < u.length(); ++ii)
          relative_mass_conservation_error = std::max(
              relative_mass_conservation_error, std::abs(initial_mass - compute_mass(u[ii].vector())) / initial_mass);
        data["quantity"]["rel. mass conserv. error"] = relative_mass_conservation_error;
      } else
        DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
                   "I do not know how to compute the quantity '" << id << "'!");
    }
    return data;
  } // ... compute(...)

protected:
  const F& flux() const override final
  {
    return Problem::access().flux;
  }

  DF make_initial_values(const S& space) override final
  {
    return Problem::access().template make_initial_values<V>(space);
  }

  GP make_initial_grid() override final
  {
    return Problem::access().make_initial_grid();
  }

  std::unique_ptr<O> make_lhs_operator(const S& space) override final
  {
    if (boundary_treatment.empty())
      return BaseType::make_lhs_operator(space);
    auto& self = *this;
    const auto& euler_tools = Problem::access().euler_tools;
    DUNE_THROW_IF(
        !self.numerical_flux_, XT::Common::Exceptions::you_are_using_this_wrong, "call set_numerical_flux() first!");
    // the layer is periodic and the operator includes handling of periodic boundaries, so we need to make an exception
    // for all non-periodic boundaries
    if (boundary_treatment == "impermeable_walls_by_direct_euler_treatment") {
      boundary_info = std::make_unique<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>>();
      boundary_info->register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
      boundary_info->register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
      XT::Grid::ApplyOn::CustomBoundaryIntersections<GV> impermeable_wall_filter(*boundary_info,
                                                                                 new XT::Grid::ImpermeableBoundary());
      auto op = std::make_unique<AdvectionFvOperator<GV, V, m>>(space.grid_view(),
                                                                *self.numerical_flux_,
                                                                space,
                                                                space,
                                                                /*periodicity_restriction=*/impermeable_wall_filter);
      // the actual handling of impermeable walls
      op->append([&](const auto& u,
                     const auto& n,
                     const auto& /*param*/ = {}) { return euler_tools.flux_at_impermeable_walls(u, n); },
                 {},
                 impermeable_wall_filter);
      return op;
    } else
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "boundary_treatment = " << boundary_treatment);
    return nullptr;
  } // ... make_lhs_operator(...)

  XT::LA::ListVectorArray<V> solve(const S& space, const double T_end, const double dt) override final
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    return solve_instationary_system_explicit_euler(u_0, *op, T_end, dt);
  }

  size_t visualization_steps_;
  std::string boundary_treatment;
  std::unique_ptr<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>> boundary_info;
}; // class InviscidCompressibleFlowEulerExplicitFvTest


using InviscidCompressibleFlow1dEulerExplicitFvTest =
    InviscidCompressibleFlowEulerExplicitFvTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest, periodic_boundaries)
{
  //  this->visualization_steps_ = 100; // <- something like this to visualize
  this->set_numerical_flux("vijayasundaram");
  this->run();
}
TEST_F(InviscidCompressibleFlow1dEulerExplicitFvTest, impermeable_walls_by_direct_euler_treatment)
{
  this->set_numerical_flux("vijayasundaram");
  this->boundary_treatment = "impermeable_walls_by_direct_euler_treatment";
  this->run();
}


#if 0
// TEST_F(InviscidCompressibleFlowEuler1dExplicitFV, impermeable_walls_by_inviscid_mirror_treatment)
//{
//  this->impermeable_walls_by_inviscid_mirror_treatment();
//}
// TEST_F(InviscidCompressibleFlowEuler1dExplicitFV,
//       inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right)
//{
//  this->inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right(10*100);
//}
#endif // 0
