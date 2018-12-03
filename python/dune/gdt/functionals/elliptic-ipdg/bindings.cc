// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#  include <dune/common/parallel/mpihelper.hh>

#  include <dune/pybindxi/pybind11.h>
#  include <dune/pybindxi/stl.h>

#  include <python/dune/xt/common/bindings.hh>
#  include <python/dune/gdt/shared.hh>

#  include <python/dune/gdt/functionals/elliptic-ipdg/bindings.hh>


PYBIND11_MODULE(__functionals_elliptic_ipdg, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  Dune::XT::Common::bindings::addbind_exceptions(m);

  py::module::import("dune.xt.common");
  py::module::import("dune.xt.grid");
  py::module::import("dune.xt.functions");
  py::module::import("dune.xt.la");
  py::module::import("dune.gdt.__spaces");

// alu_istl.cc
#  if HAVE_DUNE_ALUGRID && HAVE_DUNE_ISTL
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALU(m, leaf, view, dg, gdt, 1, istl_dense);
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALU(m, level, view, dg, gdt, 1, istl_dense);
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_ALU(m, dd_subdomain, view, dg, gdt, 1, istl_dense);
#  endif

// yasp_istl.cc
#  if HAVE_DUNE_ISTL
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_YASP(m, leaf, view, dg, gdt, 1, istl_dense);
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_YASP(m, level, view, dg, gdt, 1, istl_dense);
  DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_YASP(m, dd_subdomain, view, dg, gdt, 1, istl_dense);
#  endif

  add_initialization(m, "dune.gdt.functionals.elliptic-ipdg");
}

#endif // HAVE_DUNE_PYBINDXI
