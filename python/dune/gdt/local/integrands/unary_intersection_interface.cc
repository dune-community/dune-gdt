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

#include <dune/gdt/local/integrands/interfaces.hh>

#include <python/dune/xt/common/configuration.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/grid/grids.bindings.hh>

#include "unary_intersection_combined.hh"
#include "unary_intersection_interface.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct LocalUnaryIntersectionIntegrandInterface_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  using I = Dune::XT::Grid::extract_intersection_t<GV>;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::LocalUnaryIntersectionIntegrandInterface;
    using Dune::GDT::bindings::LocalUnaryIntersectionIntegrandSum;

    LocalUnaryIntersectionIntegrandInterface<G, I>::bind(m);
    if (d > 1) {
      LocalUnaryIntersectionIntegrandInterface<G, I, d, 1>::bind(m);
      LocalUnaryIntersectionIntegrandInterface<G, I, d, d>::bind(m);
    }
    // add your extra dimensions here
    // ...

    // need to bind LocalUnaryIntersectionIntegrandSum here due to circular dep with interface
    LocalUnaryIntersectionIntegrandSum<G, I>::bind(m);
    if (d > 1) {
      LocalUnaryIntersectionIntegrandSum<G, I, d, 1>::bind(m);
      LocalUnaryIntersectionIntegrandSum<G, I, d, d>::bind(m);
    }
    // add your extra dimensions here
    // ...

    LocalUnaryIntersectionIntegrandInterface_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct LocalUnaryIntersectionIntegrandInterface_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_local_integrands_unary_intersection_interface, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  LocalUnaryIntersectionIntegrandInterface_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
