// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/gdt/shared.hh>

#include <python/dune/gdt/functionals/elliptic-ipdg/bindings.hh>


PYBIND11_MODULE(__functionals_elliptic_ipdg, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(
      m, GDT_BINDINGS_GRID::dimension, GDT_BINDINGS_GRID, leaf, view, dg, gdt, 1, istl_dense);
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(
      m, GDT_BINDINGS_GRID::dimension, GDT_BINDINGS_GRID, level, view, dg, gdt, 1, istl_dense);
  _DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_METHODS_D(
      m, GDT_BINDINGS_GRID::dimension, GDT_BINDINGS_GRID, dd_subdomain, view, dg, gdt, 1, istl_dense);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.functionals.elliptic-ipdg");
}
