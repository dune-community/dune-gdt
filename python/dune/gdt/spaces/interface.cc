// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/xt/grid/grids.hh>

#include "interface.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct SpaceInterface_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::SpaceInterface;
    using Dune::XT::Grid::bindings::grid_name;

    SpaceInterface<GV>::bind(m, grid_name<G>::value());
    if (d > 1)
      SpaceInterface<GV, d>::bind(m, grid_name<G>::value());
    // add your extra dimensions here
    // ...
    SpaceInterface_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct SpaceInterface_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_spaces_interface, m)
{
  namespace py = pybind11;
  using namespace Dune;
  using namespace Dune::XT;
  using namespace Dune::GDT;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");

  SpaceInterface_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
