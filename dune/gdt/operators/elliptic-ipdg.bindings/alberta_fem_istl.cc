// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/gdt/operators/elliptic-ipdg.bindings.hh>


PYBIND11_PLUGIN(__operators_elliptic_ipdg_alberta_fem_istl)
{
  namespace py = pybind11;

  py::module m("__operators_elliptic_ipdg_alberta_fem_istl", "dune-gdt: EllipticMatrixOperator");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

  //  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_FEM_ISTL(m);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
