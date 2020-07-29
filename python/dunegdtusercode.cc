// This file is part of the dune-xt project:
//   https://github.com/dune-community/dune-xt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>
#include <python/dune/xt/common/parameter.hh>

#include <dune/xt/grid/gridprovider/provider.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/grid-function.hh>

using namespace Dune::XT;

PYBIND11_MODULE(dunegdtusercode, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt");

  m.def("exp_generic", [](const Grid::GridProvider<ALU_2D_SIMPLEX_CONFORMING>&, const int order) {
    return new Functions::GridFunction<Grid::extract_entity_t<ALU_2D_SIMPLEX_CONFORMING>>(
        {order, [](const auto& x, const auto&) { return std::exp(x[0] * x[1]); }});
  });
  m.def("exp_generic", [](const Grid::GridProvider<YASP_2D_EQUIDISTANT_OFFSET>&, const int order) {
    return new Functions::GridFunction<Grid::extract_entity_t<YASP_2D_EQUIDISTANT_OFFSET>>(
        {order, [](const auto& x, const auto&) { return std::exp(x[0] * x[1]); }});
  });

  // your code goes here ...
  m.def(
      "example_function", [](Common::Parameter param) { std::cout << "param = " << param << std::endl; }, "param"_a);
} // PYBIND11_MODULE(dunegdtusercode, ...)
