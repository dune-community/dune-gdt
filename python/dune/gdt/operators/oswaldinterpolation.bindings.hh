// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
#define DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>
#include <dune/xt/grid/gridprovider/provider.hh>

#include <dune/gdt/spaces.hh>

#include <dune/gdt/operators/oswaldinterpolation.hh>

namespace Dune {
namespace GDT {
namespace bindings {


template <class G,
          SpaceType space_type,
          Backends space_backend,
          XT::Grid::Layers space_layer_type,
          int p,
          class R,
          size_t r,
          XT::LA::Backends la_backend,
          XT::Grid::Layers interpolation_layer_type = space_layer_type,
          XT::Grid::Backends interpolation_layer_backend =
              SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::layer_backend>
class OswaldInterpolationOperator
{
  typedef typename SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::type S;
  typedef typename S::GridLayerType GL;
  typedef typename XT::LA::Container<R, la_backend>::VectorType V;
  typedef typename XT::Grid::
      Layer<G, interpolation_layer_type, interpolation_layer_backend, XT::Grid::DD::SubdomainGrid<G>>::type
          InterpolationLayerType;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("apply_oswald_interpolation_operator",
          [](const XT::Grid::GridProvider<G, XT::Grid::DD::SubdomainGrid<G>>& dd_grid_provider,
             const ssize_t layer_level_or_subdomain,
             const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<InterpolationLayerType>>& boundary_info,
             const GDT::ConstDiscreteFunction<S, V>& source,
             GDT::DiscreteFunction<S, V>& range) {
            GDT::OswaldInterpolationOperator<InterpolationLayerType, R>(
                dd_grid_provider.template layer<interpolation_layer_type, interpolation_layer_backend>(
                    layer_level_or_subdomain),
                boundary_info)
                .apply(source, range);
          },
          "dd_grid_provider"_a,
          "layer_level_or_subdomain"_a = -1,
          "boundary_info"_a =
              XT::Grid::AllDirichletBoundaryInfo<XT::Grid::extract_intersection_t<InterpolationLayerType>>(),
          "source"_a,
          "range"_a);
  }
}; // class OswaldInterpolationOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
