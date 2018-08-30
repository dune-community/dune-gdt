// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH
#define DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH

#include <dune/xt/common/bisect.hh>
#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/functions/lambda/function.hh>

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
    : public InstationaryEocStudy<XT::Grid::PeriodicGridLayer<typename XT::Grid::Layer<G,
                                                                                       XT::Grid::Layers::leaf,
                                                                                       XT::Grid::Backends::view>::type>,
                                  m,
                                  la>
{
  using BaseType =
      InstationaryEocStudy<XT::Grid::PeriodicGridLayer<
                               typename XT::Grid::Layer<G, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::type>,
                           m,
                           la>;

protected:
  using typename BaseType::GV;
  using typename BaseType::O;
  using typename BaseType::M;
  using typename BaseType::V;
  using typename BaseType::R;
  using typename BaseType::S;
  using typename BaseType::DF;
  using typename BaseType::GP;
  using BaseType::d;

  using F = XT::Functions::FunctionInterface<m, d, m>;
  using NF = NumericalFluxInterface<d, m>;

public:
  InstationaryNonconformingHyperbolicEocStudy(
      const double T_end,
      std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualizer =
          [](const auto& /*solution*/, const auto& /*prefix*/) { /*no visualization by default*/ })
    : BaseType(T_end, visualizer)
    , space_type_("")
    , numerical_flux_type_("")
  {
  }

  virtual std::vector<std::string> quantities() const override
  {
    auto qq = BaseType::quantities();
    if (this->space_type_ != "fv")
      qq.push_back("           CFL  "); // <- This is on purpose, see column header formatting in ConvergenceStudy.
    return qq;
  }

  virtual std::map<std::string, std::map<std::string, double>>
  compute(const size_t refinement_level,
          const std::vector<std::string>& actual_norms,
          const std::vector<std::pair<std::string, std::string>>& actual_estimates,
          const std::vector<std::string>& actual_quantities) override
  {
    auto& self = *this;
    auto quantities_to_compute = actual_quantities;
    if (self.space_type_ != "fv") {
      const auto search_result =
          std::find(quantities_to_compute.begin(), quantities_to_compute.end(), "           CFL  ");
      if (search_result != quantities_to_compute.end()) {
        self.current_data_["quantity"]["           CFL  "] =
            self.current_data_["target"]["dt"] / self.current_data_["info"]["explicit_dt"];
        quantities_to_compute.erase(search_result);
      }
    }
    auto data = BaseType::compute(refinement_level, actual_norms, actual_estimates, quantities_to_compute);
    data["quantity"]["           CFL  "] = self.current_data_["quantity"]["           CFL  "];
    return data;
  } // ... compute(...)

protected:
  virtual const F& flux() const = 0;

  virtual DF make_initial_values(const S& space) = 0;

  virtual std::pair<double, double> estimate_dt(const S& space) override
  {
    const auto u_0 = this->make_initial_values(space);
    const auto fv_dt = estimate_dt_for_hyperbolic_system(space.grid_view(), u_0, flux());
    if (space_type_ == "fv")
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
  } // ... estimate_dt(...)

  virtual std::unique_ptr<S> make_space(const GP& current_grid) override
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

  virtual std::unique_ptr<O> make_lhs_operator(const S& space) override
  {
    std::unique_ptr<NF> numerical_flux;
    if (numerical_flux_type_ == "upwind")
      numerical_flux = std::make_unique<NumericalUpwindFlux<d, m>>(flux());
    else if (numerical_flux_type_ == "vijayasundaram")
      numerical_flux = std::make_unique<NumericalVijayasundaramFlux<d, m>>(flux());
    else if (numerical_flux_type_ == "lax_friedrichs")
      numerical_flux = std::make_unique<NumericalLaxFriedrichsFlux<d, m>>(flux());
    else if (numerical_flux_type_ == "engquist_osher")
      numerical_flux = std::make_unique<NumericalEngquistOsherFlux<d, m>>(flux());
    else {
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "numerical_flux_type_ = " << numerical_flux_type_);
      return nullptr;
    }
    if (space_type_ == "fv")
      return std::make_unique<AdvectionFvOperator<M, GV, m>>(space.grid_view(), *numerical_flux, space, space);
    else
      return std::make_unique<AdvectionDgArtificialViscosityOperator<M, GV, m>>(
          space.grid_view(),
          *numerical_flux,
          space,
          space,
          /*periodicity_exception=*/XT::Grid::ApplyOn::NoIntersections<GV>(),
          DXTC_CONFIG_GET("nu_1", 0.2),
          DXTC_CONFIG_GET("alpha_1", 1.0));
  } // ... make_lhs_operator(...)

  std::string space_type_;
  std::string numerical_flux_type_;
}; // struct InstationaryNonconformingHyperbolicEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH
