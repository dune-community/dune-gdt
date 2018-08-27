// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH
#define DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH

#include <cmath>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/functions/base/sliced.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/local/integrands/identity.hh>
#include <dune/gdt/test/instationary-eocstudies/hyperbolic-nonconforming.hh>
#include <dune/gdt/tools/euler.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct InviscidCompressibleFlowEulerProblem
{
  static const constexpr size_t d = G::dimension;
  static const constexpr size_t m = d + 2;
  using DomainType = XT::Common::FieldVector<double, d>;
  using RangeType = XT::Common::FieldVector<double, m>;

  const EulerTools<d> euler_tools;
  const XT::Functions::LambdaFunction<m, d, m> flux;
  const double T_end;

  InviscidCompressibleFlowEulerProblem()
    : euler_tools(1.4) // air or water at roughly 20 deg Cels.
    , flux(euler_tools.flux_order(),
           [&](const auto& u, const auto& /*param*/) { return euler_tools.flux(u); },
           "euler_flux",
           {},
           [&](const auto& u, const auto& /*param*/) { return euler_tools.flux_jacobian(u); })
    , T_end(1.)
  {
  }

  XT::Grid::GridProvider<G> make_initial_grid() const
  {
    return XT::Grid::make_cube_grid<G>(-1., 1., 16u);
  }

  template <class Vector, class GV>
  DiscreteFunction<Vector, GV, m> make_initial_values(const SpaceInterface<GV, m>& space) const
  {
    return interpolate<Vector>(
        0,
        [&](const auto& xx, const auto& /*mu*/) {
          RangeType primitive_variables(0.);
          // density
          if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
            primitive_variables[0] = 4.;
          else
            primitive_variables[0] = 1.;
          // velocity
          for (size_t ii = 0; ii < d; ++ii)
            primitive_variables[1 + ii] = 0.;
          // pressure
          if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
            primitive_variables[m - 1] = 1.6;
          else
            primitive_variables[m - 1] = 0.4;
          return euler_tools.to_conservative(primitive_variables);
        },
        space);
  } // ... make_initial_values(...)
}; // struct InviscidCompressibleFlowEulerProblem


template <class G>
class InviscidCompressibleFlowEulerTest
    : protected XT::Common::ConstStorageProvider<InviscidCompressibleFlowEulerProblem<G>>,
      public InstationaryNonconformingHyperbolicEocStudy<G, G::dimension + 2>
{
  using BaseType = InstationaryNonconformingHyperbolicEocStudy<G, G::dimension + 2>;

protected:
  using Problem = XT::Common::ConstStorageProvider<InviscidCompressibleFlowEulerProblem<G>>;
  using BaseType::d;
  using BaseType::m;
  using typename BaseType::E;
  using typename BaseType::F;
  using typename BaseType::R;
  using typename BaseType::DF;
  using typename BaseType::GV;
  using typename BaseType::GP;
  using typename BaseType::S;
  using typename BaseType::BS;
  using typename BaseType::V;
  using typename BaseType::O;

public:
  InviscidCompressibleFlowEulerTest()
    : Problem(new InviscidCompressibleFlowEulerProblem<G>())
    , BaseType(Problem::access().T_end,
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

protected:
  const F& flux() const override final
  {
    return Problem::access().flux;
  }

  DF make_initial_values(const S& space) override final
  {
    if (boundary_treatment == "inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right") {
      const auto& euler_tools = Problem::access().euler_tools;
      return interpolate<V>(0,
                            [&](const auto& /*xx*/, const auto& /*param*/) {
                              FieldVector<R, m> primitive_variables(0.);
                              // density
                              primitive_variables[0] = 0.5;
                              // velocity
                              for (size_t ii = 0; ii < d; ++ii)
                                primitive_variables[1 + ii] = 0.;
                              // pressure
                              primitive_variables[m - 1] = 0.4;
                              return euler_tools.to_conservative(primitive_variables);
                            },
                            space);
    } else
      return Problem::access().template make_initial_values<V>(space);
  } // ... make_initial_values(...)

  // Apart from the boundary_data_range, this is the same as BaseType::estimate_dt
  std::pair<double, double> estimate_dt(const S& space) override final
  {
    if (boundary_treatment == "inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right") {
      const auto u_0 = this->make_initial_values(space);
      const auto fv_dt = estimate_dt_for_hyperbolic_system(space.grid_view(),
                                                           u_0,
                                                           flux(),
                                                           /*boundary_data_range=*/{{0.5, 0., 0.4}, {1.5, 0.5, 0.4}});
      if (this->space_type_ == "fv")
        return {fv_dt, fv_dt};
      const auto max_sup_norm = 1.25 * u_0.dofs().vector().sup_norm();
      const auto actual_dt = XT::Common::find_largest_by_bisection(1e-15, fv_dt, [&](const auto& dt_to_test) {
        const auto solution = this->solve(space, 250 * dt_to_test, dt_to_test);
        for (const auto& vec : solution.vectors())
          if (vec.sup_norm() > max_sup_norm)
            return false;
        return true;
      });
      return {fv_dt, actual_dt};
    } else
      return BaseType::estimate_dt(space);
  } // ... estimate_dt(...)

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
    auto numerical_flux = make_numerical_vijayasundaram_flux(self.flux());
    // the layer is periodic and the operator includes handling of periodic boundaries, so we need to make an exception
    // for all non-periodic boundaries
    if (boundary_treatment == "impermeable_walls_by_direct_euler_treatment") {
      boundary_info = std::make_unique<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>>();
      boundary_info->register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
      boundary_info->register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
      const XT::Grid::ApplyOn::CustomBoundaryIntersections<GV> impermeable_wall_filter(
          *boundary_info, new XT::Grid::ImpermeableBoundary());
      auto op = std::make_unique<AdvectionFvOperator<GV, V, m>>(space.grid_view(),
                                                                numerical_flux,
                                                                space,
                                                                space,
                                                                /*periodicity_restriction=*/impermeable_wall_filter);
      // the actual handling of impermeable walls
      op->append(/*numerical_boundary_flux=*/[&](
                     const auto& u,
                     const auto& n,
                     const auto& /*param*/) { return euler_tools.flux_at_impermeable_walls(u, n); },
                 {},
                 impermeable_wall_filter);
      return op;
    } else if (boundary_treatment == "impermeable_walls_by_inviscid_mirror_treatment") {
      boundary_info = std::make_unique<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>>();
      boundary_info->register_new_normal({-1}, new XT::Grid::ImpermeableBoundary());
      boundary_info->register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
      const XT::Grid::ApplyOn::CustomBoundaryIntersections<GV> impermeable_wall_filter(
          *boundary_info, new XT::Grid::ImpermeableBoundary());
      auto op = std::make_unique<AdvectionFvOperator<GV, V, m>>(space.grid_view(),
                                                                numerical_flux,
                                                                space,
                                                                space,
                                                                /*periodicity_restriction=*/impermeable_wall_filter);
      // the actual handling of impermeable walls, see [DF2015, p. 415, (8.66 - 8.67)]
      op->append(
          /*boundary_extrapolation=*/
          [&](const auto& intersection,
              const auto& xx_in_reference_intersection_coordinates,
              const auto& /*flux*/,
              const auto& u,
              const auto& /*param*/) {
            const auto normal = intersection.unitOuterNormal(xx_in_reference_intersection_coordinates);
            const auto rho = euler_tools.density_from_conservative(u);
            auto velocity = euler_tools.velocity_from_conservative(u);
            velocity -= normal * 2. * (velocity * normal);
            const auto pressure = euler_tools.pressure_from_conservative(u);
            return euler_tools.to_conservative(XT::Common::hstack(rho, velocity, pressure));
          },
          {},
          impermeable_wall_filter);
      return op;
    } else if (boundary_treatment == "inflow_from_the_left_by_heuristic_euler_treatment_impermeable_wall_right") {
      periodic_density_variation_ =
          std::unique_ptr<XT::Functions::LambdaFunction<d, m>>(new XT::Functions::LambdaFunction<d, m>(
              /*order=*/0,
              [&](const auto& /*xx*/, const auto& param) {
                FieldVector<R, m> primitive_variables(0.);
                const auto t = param.get("_t").at(0);
                // density
                primitive_variables[0] = 1. + 0.5 * std::sin(2. * M_PI * (t - 0.25));
                // velocity
                for (size_t ii = 0; ii < d; ++ii)
                  primitive_variables[1 + ii] = 0.25 + 0.25 * std::sin(2. * M_PI * (t - 0.25));
                // pressure
                primitive_variables[m - 1] = 0.4;
                return euler_tools.to_conservative(primitive_variables);
              },
              /*name=*/"periodic_density_variation",
              /*parameter_type=*/{"_t", 1}));
      local_periodic_density_variation_ = periodic_density_variation_->template as_grid_function<E>().local_function();
      // inflow/outflow left, impermeable wall right
      boundary_info = std::make_unique<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>>();
      boundary_info->register_new_normal({-1}, new XT::Grid::InflowOutflowBoundary());
      boundary_info->register_new_normal({1}, new XT::Grid::ImpermeableBoundary());
      const XT::Grid::ApplyOn::CustomBoundaryIntersections<GV> impermeable_wall_filter(
          *boundary_info, new XT::Grid::ImpermeableBoundary());
      const XT::Grid::ApplyOn::CustomBoundaryIntersections<GV> inflow_outflow_filter(
          *boundary_info, new XT::Grid::InflowOutflowBoundary());
      auto op = std::make_unique<AdvectionFvOperator<GV, V, m>>(
          space.grid_view(),
          numerical_flux,
          space,
          space,
          /*periodicity_restriction=*/*(inflow_outflow_filter || impermeable_wall_filter));
      // the actual handling of inflow/outflow, see [DF2015, p. 421, (8.88)]
      // (supposedly this does not work well for slow flows!)
      const auto heuristic_euler_inflow_outflow_treatment = [&](const auto& intersection,
                                                                const auto& xx_in_reference_intersection_coordinates,
                                                                const auto& /*flux*/,
                                                                const auto& u,
                                                                const auto& param) {
        using RangeType = XT::Common::FieldVector<R, m>;
        // evaluate boundary values
        const auto element = intersection.inside();
        const auto xx_in_reference_element_coordinates =
            intersection.geometryInInside().global(xx_in_reference_intersection_coordinates);
        local_periodic_density_variation_->bind(element);
        const RangeType bv = local_periodic_density_variation_->evaluate(xx_in_reference_element_coordinates, param);
        // determine flow regime
        const auto a = euler_tools.speed_of_sound_from_conservative(u);
        const auto velocity = euler_tools.velocity_from_conservative(u);
        const auto normal = intersection.unitOuterNormal(xx_in_reference_intersection_coordinates);
        const auto flow_speed = velocity * normal;
        // compute v
        if (flow_speed < -a) {
          // supersonic inlet
          return bv;
        } else if (!(flow_speed > 0)) {
          // subsonic inlet
          const auto rho_outer = euler_tools.density_from_conservative(bv);
          const auto v_outer = euler_tools.velocity_from_conservative(bv);
          const auto p_inner = euler_tools.pressure_from_conservative(u);
          return euler_tools.to_conservative(XT::Common::hstack(rho_outer, v_outer, p_inner));
        } else if (flow_speed < a) {
          // subsonic outlet
          const auto rho_inner = euler_tools.density_from_conservative(u);
          const auto v_inner = euler_tools.velocity_from_conservative(u);
          const auto p_outer = euler_tools.pressure_from_conservative(bv);
          return euler_tools.to_conservative(XT::Common::hstack(rho_inner, v_inner, p_outer));
        } else {
          // supersonic outlet
          return RangeType(u);
        }
      }; // ... heuristic_euler_inflow_outflow_treatment(...)
      op->append(heuristic_euler_inflow_outflow_treatment, {}, inflow_outflow_filter);
      // the actual handling of impermeable walls, see above
      op->append(
          /*boundary_extrapolation=*/
          [&](const auto& intersection,
              const auto& xx_in_reference_intersection_coordinates,
              const auto& /*flux*/,
              const auto& u,
              const auto& /*param*/) {
            const auto normal = intersection.unitOuterNormal(xx_in_reference_intersection_coordinates);
            const auto rho = euler_tools.density_from_conservative(u);
            auto velocity = euler_tools.velocity_from_conservative(u);
            velocity -= normal * 2. * (velocity * normal);
            const auto pressure = euler_tools.pressure_from_conservative(u);
            return euler_tools.to_conservative(XT::Common::hstack(rho, velocity, pressure));
          },
          {},
          impermeable_wall_filter);
      return op;
    } else
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "boundary_treatment = " << boundary_treatment);
    return nullptr;
  } // ... make_lhs_operator(...)

  size_t visualization_steps_;
  std::string boundary_treatment;
  std::unique_ptr<XT::Grid::NormalBasedBoundaryInfo<XT::Grid::extract_intersection_t<GV>>> boundary_info;
  std::unique_ptr<XT::Functions::LambdaFunction<d, m>> periodic_density_variation_;
  std::unique_ptr<typename XT::Functions::ElementFunctionInterface<E, m>> local_periodic_density_variation_;
}; // class InviscidCompressibleFlowEulerTest


template <class G>
class InviscidCompressibleFlowEulerExplicitTest : public InviscidCompressibleFlowEulerTest<G>
{
  using BaseType = InviscidCompressibleFlowEulerTest<G>;
  using typename BaseType::S;
  using typename BaseType::V;

protected:
  XT::LA::ListVectorArray<V> solve(const S& space, const double T_end, const double dt) override final
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    return solve_instationary_system_explicit_euler(u_0, *op, T_end, dt);
  }
}; // class InviscidCompressibleFlowEulerExplicitTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH
