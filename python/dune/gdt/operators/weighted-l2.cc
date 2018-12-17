// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   Rene Milk       (2018)
//   Tobias Leibner  (2017)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/available_types.hh>
#include <python/dune/gdt/shared.hh>
#include <python/dune/xt/grid/grids.bindings.hh>
#include <dune/xt/grid/layers.hh>

#include <python/dune/gdt/operators/weighted-l2.hh>


PYBIND11_MODULE(__operators_weighted_l2, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;
  using Dune::XT::Grid::Backends;
  using Dune::XT::Grid::Layers;

  Dune::GDT::bindings::WeightedL2LocalizableProduct<GDT_BINDINGS_GRID, Layers::leaf, Backends::view>::bind(m);
  Dune::GDT::bindings::WeightedL2LocalizableProduct<GDT_BINDINGS_GRID, Layers::level, Backends::view>::bind(m);
  Dune::GDT::bindings::WeightedL2LocalizableProduct<GDT_BINDINGS_GRID, Layers::dd_subdomain, Backends::view>::bind(m);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.operators.elliptic");
}
