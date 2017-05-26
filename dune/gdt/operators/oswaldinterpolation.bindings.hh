// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
#define DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include "oswaldinterpolation.hh"

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
          XT::LA::Backends la_backend /*,
          XT::Grid::Layers layer_type = space_layer_type,
          XT::Grid::Backends layer_backend =
              SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::layer_backend*/>
class OswaldInterpolationOperator
{
  typedef typename SpaceProvider<G, space_layer_type, space_type, space_backend, p, R, r>::type S;
  typedef typename S::GridLayerType GL;
  typedef typename XT::LA::Container<R, la_backend>::VectorType V;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("apply_oswald_interpolation_operator",
          [](const GDT::ConstDiscreteFunction<S, V>& source,
             GDT::DiscreteFunction<S, V>& range,
             const bool zero_boundary) {
            GDT::OswaldInterpolationOperator<GL, R>(source.space().grid_layer(), zero_boundary).apply(source, range);
          },
          "source"_a,
          "range"_a,
          "zero_boundary"_a = true);
  }
}; // class OswaldInterpolationOperator


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_OPERATORS_OSWALDINTERPOLATION_BINDINGS_HH
