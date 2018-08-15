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
      const size_t num_refinements = 3,
      const size_t num_additional_refinements_for_reference = 2,
      std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualizer =
          [](const auto& /*solution*/, const auto& /*prefix*/) { /*no visualization by default*/ })
    : BaseType(T_end, num_refinements, num_additional_refinements_for_reference, visualizer)
    , space_type_("")
    , numerical_flux_type_("")
  {
  }

protected:
  virtual const F& flux() const = 0;

  virtual DF make_initial_values(const S& space) = 0;

  virtual double estimate_dt(const S& space) override
  {
    return estimate_dt_for_hyperbolic_system(space.grid_view(), make_initial_values(space), flux());
  }

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
      return std::make_unique<AdvectionFvOperator<GV, V, m>>(space.grid_view(), *numerical_flux, space, space);
    else
      return std::make_unique<AdvectionDgOperator<GV, V, m>>(space.grid_view(), *numerical_flux, space, space);
  } // ... make_lhs_operator(...)

  std::string space_type_;
  std::string numerical_flux_type_;
}; // struct InstationaryNonconformingHyperbolicEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_NONCONFORMING_HH
