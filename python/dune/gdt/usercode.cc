// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#include "config.h"

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>
#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/common/fvector.hh>
#include <python/dune/xt/common/exceptions.bindings.hh>


PYBIND11_MODULE(usercode, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.la");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.gdt.discretefunction");

  Dune::XT::Common::bindings::addbind_exceptions(m);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt");

  // put your bindings below
} // PYBIND11_MODULE(usercode, ...)
