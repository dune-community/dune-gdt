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

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <dune/xt/common/bindings.hh>

#include <dune/gdt/prolongations.bindings.hh>


PYBIND11_PLUGIN(__prolongations)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__prolongations", "dune-gdt: prolongations");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.prolongations");


  DUNE_GDT_PROLONGATIONS_BIND(m);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
