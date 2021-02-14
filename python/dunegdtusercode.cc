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


PYBIND11_MODULE(dunegdtusercode, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt");

  // your code goes here ...
  m.def(
      "example_function",
      [](Dune::XT::Common::Parameter param) { std::cout << "param = " << param << std::endl; },
      "param"_a);
} // PYBIND11_MODULE(dunegdtusercode, ...)
