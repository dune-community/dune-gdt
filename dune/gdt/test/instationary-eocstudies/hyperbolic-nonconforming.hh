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

#ifndef DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH
#define DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH

#include <dune/xt/common/bisect.hh>
#include <dune/xt/test/common.hh>

#include <dune/xt/grid/view/periodic.hh>

#include <dune/xt/functions/generic/function.hh>

#include <dune/gdt/local/numerical-fluxes/engquist-osher.hh>
#include <dune/gdt/local/numerical-fluxes/lax-friedrichs.hh>
#include <dune/gdt/local/numerical-fluxes/upwind.hh>
#include <dune/gdt/local/numerical-fluxes/vijayasundaram.hh>
#include <dune/gdt/operators/advection-dg.hh>
#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/hyperbolic.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class G, size_t m = 1, XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class InstationaryNonconformingHyperbolicEocStudy
  : public InstationaryEocStudy<
        XT::Grid::PeriodicGridLayer<
            typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>,
        m,
        la>
{
  using BaseType = InstationaryEocStudy<
      XT::Grid::PeriodicGridLayer<typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>,
      m,
      la>;

protected:
  static const constexpr int m_as_int = int(ssize_t(m));
  using BaseType::d;
  using typename BaseType::DF;
  using typename BaseType::GP;
  using typename BaseType::GV;
  using typename BaseType::M;
  using typename BaseType::O;
  using typename BaseType::R;
  using typename BaseType::S;
  using typename BaseType::V;
  using I = XT::Grid::extract_intersection_t<GV>;

  using F = XT::Functions::FunctionInterface<m, d, m>;
  using NF = NumericalFluxInterface<I, d, m>;

public:
  InstationaryNonconformingHyperbolicEocStudy(
      const double T_end,
      const std::string timestepping,
      std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualizer =
          [](const auto& /*solution*/, const auto& /*prefix*/) { /*no visualization by default*/ },
      const size_t num_refinements = DXTC_TEST_CONFIG_GET("setup.num_refinements", 3),
      const size_t num_additional_refinements_for_reference =
          DXTC_TEST_CONFIG_GET("setup.num_additional_refinements_for_reference", 2))
    : BaseType(T_end, timestepping, visualizer, num_refinements, num_additional_refinements_for_reference)
    , space_type_("")
    , numerical_flux_type_("")
    , dg_artificial_viscosity_nu_1_(
          DXTC_TEST_CONFIG_GET("setup.dg_artificial_viscosity_nu_1", advection_dg_artificial_viscosity_default_nu_1()))
    , dg_artificial_viscosity_alpha_1_(DXTC_TEST_CONFIG_GET("setup.dg_artificial_viscosity_alpha_1",
                                                            advection_dg_artificial_viscosity_default_alpha_1()))
    , dg_artificial_viscosity_component_(DXTC_TEST_CONFIG_GET("setup.dg_artificial_viscosity_component",
                                                              advection_dg_artificial_viscosity_default_component()))
  {}

protected:
  virtual const F& flux() const = 0;

  virtual DF make_initial_values(const S& space) = 0;

  std::unique_ptr<S> make_space(const GP& current_grid) override
  {
    if (space_type_ == "fv")
      return std::make_unique<FiniteVolumeSpace<GV, m>>(XT::Grid::make_periodic_grid_layer(current_grid.leaf_view()));
    else if (space_type_.size() >= 4 && space_type_.substr(0, 4) == "dg_p") {
      const auto order = XT::Common::from_string<int>(space_type_.substr(4));
      return std::make_unique<DiscontinuousLagrangeSpace<GV, m>>(
          XT::Grid::make_periodic_grid_layer(current_grid.leaf_view()), order);
    } else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "space_type_ = " << space_type_);
      return nullptr;
    }
  } // ... make_space(...)

  std::unique_ptr<O> make_lhs_operator(const S& space) override
  {
    std::unique_ptr<NF> numerical_flux;
    if (numerical_flux_type_ == "upwind")
      numerical_flux = std::make_unique<NumericalUpwindFlux<I, d, m>>(flux());
    else if (numerical_flux_type_ == "vijayasundaram")
      numerical_flux = std::make_unique<NumericalVijayasundaramFlux<I, d, m>>(flux());
    else if (numerical_flux_type_ == "lax_friedrichs")
      numerical_flux = std::make_unique<NumericalLaxFriedrichsFlux<I, d, m>>(flux());
    else if (numerical_flux_type_ == "engquist_osher")
      numerical_flux = std::make_unique<NumericalEngquistOsherFlux<I, d, m>>(flux());
    else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "numerical_flux_type_ = " << numerical_flux_type_);
      return nullptr;
    }
    if (space_type_ == "fv")
      return std::make_unique<AdvectionFvOperator<M, GV, m>>(space.grid_view(), *numerical_flux, space, space);
    else
      return std::make_unique<AdvectionDgOperator<M, GV, m>>(
          space.grid_view(),
          *numerical_flux,
          space,
          space,
          /*periodicity_exception=*/XT::Grid::ApplyOn::NoIntersections<GV>(),
          dg_artificial_viscosity_nu_1_,
          dg_artificial_viscosity_alpha_1_,
          dg_artificial_viscosity_component_);
  } // ... make_lhs_operator(...)

  virtual double estimate_fixed_explicit_fv_dt(
      const S& space,
      const std::pair<XT::Common::FieldVector<R, m_as_int>, XT::Common::FieldVector<R, m_as_int>>& boundary_data_range =
          {XT::Common::FieldVector<R, m_as_int>(std::numeric_limits<R>::max()),
           XT::Common::FieldVector<R, m_as_int>(std::numeric_limits<R>::min())})
  {
    return estimate_dt_for_hyperbolic_system(
        space.grid_view(), make_initial_values(space), flux(), boundary_data_range);
  }

  virtual double
  estimate_fixed_explicit_dt(const S& space, const double max_overshoot = 1.25, const int max_steps_to_try = 250)
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    const auto max_sup_norm = max_overshoot * u_0.dofs().vector().sup_norm();
    return XT::Common::find_largest_by_bisection(
        /*min_dt=*/10 * std::numeric_limits<double>::epsilon(),
        /*max_dt=*/this->T_end_,
        /*success=*/
        [&](const auto& dt_to_test) {
          try {
            auto u = u_0.dofs().vector();
            const double T_end = max_steps_to_try * dt_to_test;
            double time = 0.;
            // explicit euler
            while (time < T_end + dt_to_test) {
              u -= op->apply(u, {{"_t", {time}}, {"_dt", {dt_to_test}}}) * dt_to_test;
              time += dt_to_test;
              if (u.sup_norm() > max_sup_norm)
                return false;
            }
            return true;
          } catch (...) {
            return false;
          }
        },
        1e-2);
  } // ... estimate_fixed_explicit_dt(...)

  virtual double estimate_fixed_explicit_dt_to_T_end(const S& space,
                                                     const double& min_dt,
                                                     const double& T_end,
                                                     const double max_overshoot = 1.25)
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    const auto max_sup_norm = max_overshoot * u_0.dofs().vector().sup_norm();
    return XT::Common::find_largest_by_bisection(
        /*min_dt=*/min_dt,
        /*max_dt=*/T_end,
        /*success=*/
        [&](const auto& dt_to_test) {
          try {
            auto u = u_0.dofs().vector();
            double time = 0.;
            // explicit euler
            while (time < T_end + dt_to_test) {
              u -= op->apply(u, {{"_t", {time}}, {"_dt", {dt_to_test}}}) * dt_to_test;
              time += dt_to_test;
              if (u.sup_norm() > max_sup_norm)
                return false;
            }
            return true;
          } catch (...) {
            return false;
          }
        },
        1e-2);
  } // ... estimate_fixed_explicit_dt(...)

  std::string space_type_;
  std::string numerical_flux_type_;
  double dg_artificial_viscosity_nu_1_;
  double dg_artificial_viscosity_alpha_1_;
  size_t dg_artificial_viscosity_component_;
}; // struct InstationaryNonconformingHyperbolicEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH
