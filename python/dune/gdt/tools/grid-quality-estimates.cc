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
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/tools/grid-quality-estimates.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/parameter.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <python/dune/xt/la/traits.hh>


namespace Dune {
namespace GDT {
namespace bindings {


template <class GV, size_t r = 1>
class estimate_inverse_inequality_constant
{
  using S = SpaceInterface<GV, r>;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("estimate_inverse_inequality_constant",
          [](const S& space) { return GDT::estimate_inverse_inequality_constant(space); },
          "space"_a,
          py::call_guard<py::gil_scoped_release>());
  } // ... bind(...)
}; // class estimate_inverse_inequality_constant


template <class GV, size_t r = 1>
class estimate_combined_inverse_trace_inequality_constant
{
  using S = SpaceInterface<GV, r>;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("estimate_combined_inverse_trace_inequality_constant",
          [](const S& space) { return GDT::estimate_combined_inverse_trace_inequality_constant(space); },
          "space"_a,
          py::call_guard<py::gil_scoped_release>());
  } // ... bind(...)
}; // class estimate_combined_inverse_trace_inequality_constant


template <class GV>
class estimate_element_to_intersection_equivalence_constant
{
  using G = XT::Grid::extract_grid_t<GV>;
  using GP = XT::Grid::GridProvider<G>;

public:
  static void bind(pybind11::module& m)
  {
    namespace py = pybind11;
    using namespace pybind11::literals;

    m.def("estimate_element_to_intersection_equivalence_constant",
          [](const GP& grid) { return GDT::estimate_element_to_intersection_equivalence_constant(grid.leaf_view()); },
          "grid_provider"_a,
          py::call_guard<py::gil_scoped_release>());
  } // ... bind(...)
}; // class estimate_element_to_intersection_equivalence_constant


} // namespace bindings
} // namespace GDT
} // namespace Dune


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct estimate_inverse_inequality_constant_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::estimate_inverse_inequality_constant<GV>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::estimate_inverse_inequality_constant<GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    estimate_inverse_inequality_constant_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct estimate_inverse_inequality_constant_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct estimate_combined_inverse_trace_inequality_constant_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::estimate_combined_inverse_trace_inequality_constant<GV>::bind(m);
    if (d > 1)
      Dune::GDT::bindings::estimate_combined_inverse_trace_inequality_constant<GV, d>::bind(m);
    // add your extra dimensions here
    // ...
    estimate_combined_inverse_trace_inequality_constant_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct estimate_combined_inverse_trace_inequality_constant_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct estimate_element_to_intersection_equivalence_constant_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    Dune::GDT::bindings::estimate_element_to_intersection_equivalence_constant<GV>::bind(m);

    estimate_element_to_intersection_equivalence_constant_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct estimate_element_to_intersection_equivalence_constant_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_tools_grid_quality_estimates, m)
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

  estimate_inverse_inequality_constant_for_all_grids<>::bind(m);
  estimate_combined_inverse_trace_inequality_constant_for_all_grids<>::bind(m);
  estimate_element_to_intersection_equivalence_constant_for_all_grids<>::bind(m);
}
