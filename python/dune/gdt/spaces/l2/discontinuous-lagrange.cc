// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <python/dune/xt/grid/grids.bindings.hh>

#include "discontinuous-lagrange.hh"


template <class GridTypes = Dune::XT::Grid::bindings::AvailableGridTypes>
struct DiscontinuousLagrangeSpace_for_all_grids
{
  using G = Dune::XT::Common::tuple_head_t<GridTypes>;
  using GV = typename G::LeafGridView;
  static const constexpr size_t d = G::dimension;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::DiscontinuousLagrangeSpace;
    using Dune::XT::Grid::bindings::grid_name;

    DiscontinuousLagrangeSpace<GV>::bind(m, grid_name<G>::value());
    if (d > 1)
      DiscontinuousLagrangeSpace<GV, d>::bind(m, grid_name<G>::value());
    // add your extra dimensions here
    // ...
    DiscontinuousLagrangeSpace_for_all_grids<Dune::XT::Common::tuple_tail_t<GridTypes>>::bind(m);
  }
};

template <>
struct DiscontinuousLagrangeSpace_for_all_grids<Dune::XT::Common::tuple_null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_spaces_l2_discontinuous_lagrange, m)
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

  DiscontinuousLagrangeSpace_for_all_grids<XT::Grid::bindings::AvailableGridTypes>::bind(m);
}
