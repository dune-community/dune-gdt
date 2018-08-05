// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_FV_HH
#define DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_FV_HH

#include <dune/xt/grid/view/periodic.hh>
#include <dune/xt/functions/lambda/function.hh>

#include <dune/gdt/operators/advection-fv.hh>
#include <dune/gdt/spaces/l2/finite-volume.hh>
#include <dune/gdt/tools/hyperbolic.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class G, size_t m = 1, XT::LA::Backends la = XT::LA::Backends::istl_sparse>
class InstationaryHyperbolicFiniteVolumeEocStudy
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
  InstationaryHyperbolicFiniteVolumeEocStudy(
      const double T_end,
      const size_t num_refinements = 3,
      const size_t num_additional_refinements_for_reference = 2,
      std::function<void(const DiscreteBochnerFunction<V, GV, m>&, const std::string&)> visualizer =
          [](const auto& /*solution*/, const auto& /*prefix*/) { /*no visualization by default*/ })
    : BaseType(T_end, num_refinements, num_additional_refinements_for_reference, visualizer)
  {
  }

protected:
  virtual const F& flux() const = 0;

  virtual DF make_initial_values(const S& space) = 0;

  virtual double estimate_dt(const S& space) override
  {
    return estimate_dt_for_hyperbolic_system(space.grid_view(), make_initial_values(space), flux());
  }

  void set_numerical_flux(const std::string type)
  {
    if (type == "upwind")
      numerical_flux_ = std::make_unique<NumericalUpwindFlux<d>>(flux());
    else if (type == "vijayasundaram")
      numerical_flux_ = std::make_unique<NumericalVijayasundaramFlux<d>>(flux());
    else if (type == "lax_friedrichs")
      numerical_flux_ = std::make_unique<NumericalLaxFriedrichsFlux<d>>(flux());
    else if (type == "engquist_osher")
      numerical_flux_ = std::make_unique<NumericalEngquistOsherFlux<d>>(flux());
    else
      DUNE_THROW(XT::Common::Exceptions::wrong_input_given, "type = " << type);
  } // ... set_numerical_flux(...)

  virtual std::unique_ptr<S> make_space(const GP& current_grid) override
  {
    return std::make_unique<FiniteVolumeSpace<GV, m>>(
        make_finite_volume_space<m>(XT::Grid::make_periodic_grid_layer(current_grid.leaf_view())));
  }

  virtual std::unique_ptr<O> make_lhs_operator(const S& space) override
  {
    DUNE_THROW_IF(
        !numerical_flux_, XT::Common::Exceptions::you_are_using_this_wrong, "call set_numerical_flux() first!");
    return std::make_unique<AdvectionFvOperator<GV, V, m>>(space.grid_view(), *numerical_flux_, space, space);
  }

  const std::string numerical_flux_type_;
  std::unique_ptr<NF> numerical_flux_;
}; // struct InstationaryHyperbolicFiniteVolumeEocStudy


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INSTATIONARY_EOCSTUDIES_HYPERBOLIC_FV_HH
