// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>

#include <python/dune/gdt/discretefunction/bindings.hh>
#include <python/dune/gdt/shared.hh>


PYBIND11_MODULE(__shared, m)
{
  namespace py = pybind11;
  m.attr("GDT_BINDINGS_GRID") = Dune::XT::Grid::bindings::grid_name<GDT_BINDINGS_GRID>::value();
}
