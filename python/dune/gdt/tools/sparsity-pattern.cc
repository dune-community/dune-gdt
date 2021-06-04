// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/gridprovider/coupling.hh>
#include <dune/xt/grid/dd/glued.hh>
#include <dune/xt/grid/view/coupling.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/tools/sparsity-pattern.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, size_t r = 1>
class make_sparsity_patterns
{
  using G = std::decay_t<XT::Grid::extract_grid_t<GV>>;
  static const size_t d = G::dimension;
  using S = SpaceInterface<GV, r>;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def(
        "make_element_sparsity_pattern",
        [](const S& space) { return make_element_sparsity_pattern(space); },
        "space"_a,
        py::call_guard<py::gil_scoped_release>());
    m.def(
        "make_intersection_sparsity_pattern",
        [](const S& space) { return make_intersection_sparsity_pattern(space); },
        "space"_a,
        py::call_guard<py::gil_scoped_release>());
    m.def(
        "make_element_and_intersection_sparsity_pattern",
        [](const S& space) { return make_element_and_intersection_sparsity_pattern(space); },
        "space"_a,
        py::call_guard<py::gil_scoped_release>());

#if HAVE_DUNE_GRID_GLUE
    if constexpr (d == 2) {
      using GridGlueType = Dune::XT::Grid::DD::Glued<G, G, Dune::XT::Grid::Layers::leaf>;
      using CGV = Dune::XT::Grid::CouplingGridView<GridGlueType>;
      using CGP = Dune::XT::Grid::CouplingGridProvider<CGV>;
      m.def(
          "make_coupling_sparsity_pattern",
          [](const S& space_in, const S& space_out, const CGP& coupling_grid) {
            return make_coupling_sparsity_pattern(space_in, space_out, coupling_grid.coupling_view());
          },
          "space_in"_a,
          "space_out"_a,
          "coupling_grid"_a,
          py::call_guard<py::gil_scoped_release>());
    }
#endif
  } // ... bind(...)
}; // class make_sparsity_patterns


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct make_sparsity_patterns_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::make_sparsity_patterns<GV>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::make_sparsity_patterns<GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    make_sparsity_patterns_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct make_sparsity_patterns_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_tools_sparsity_pattern, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  py::module::import("dune.gdt._spaces_interface");

  make_sparsity_patterns_for_all_grids<>::bind(m);
}
