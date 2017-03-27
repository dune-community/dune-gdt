// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/pybindxi/pybind11.h>

#include <dune/gdt/operators/elliptic-ipdg.bindings.hh>


PYBIND11_PLUGIN(operators_elliptic_ipdg_alu_fem_istl)
{
  namespace py = pybind11;

  py::module m("operators_elliptic_ipdg_alu_fem_istl", "dune-gdt: EllipticMatrixOperator");

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");

#if HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, leaf, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, level, part, cg, fem, 1, istl_sparse);
  DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_ALU(m, dd_subdomain, part, cg, fem, 1, istl_sparse);
#endif

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
