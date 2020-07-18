// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include "continuous-lagrange.hh"


template <class GridTypes = Dune::XT::Grid::AvailableGridTypes>
struct ContinuousLagrangeSpace_for_all_grids
{
  using G = typename GridTypes::head_type;
  using GV = typename G::LeafGridView;

  static void bind(pybind11::module& m)
  {
    using Dune::GDT::bindings::ContinuousLagrangeSpace;
    using Dune::XT::Grid::bindings::grid_name;

    ContinuousLagrangeSpace<GV>::bind(m, grid_name<G>::value());
    // add your extra dimensions here
    // ...
    ContinuousLagrangeSpace_for_all_grids<typename GridTypes::tail_type>::bind(m);
  }
};

template <>
struct ContinuousLagrangeSpace_for_all_grids<boost::tuples::null_type>
{
  static void bind(pybind11::module& /*m*/) {}
};


PYBIND11_MODULE(_spaces_h1_continuous_lagrange, m)
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

  ContinuousLagrangeSpace_for_all_grids<XT::Grid::AvailableGridTypes>::bind(m);
}
