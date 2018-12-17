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
#include <python/dune/gdt/shared.hh>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <python/dune/xt/grid/available_types.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/spaces.hh>
#include <python/dune/gdt/operators/fluxreconstruction.hh>


PYBIND11_MODULE(__operators_fluxreconstruction, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;


  Dune::GDT::bindings::DiffusiveFluxReconstructionOperator<GDT_BINDINGS_GRID,
                                                           Dune::GDT::SpaceType::rt,
                                                           Dune::GDT::Backends::gdt,
                                                           Dune::XT::Grid::Layers::leaf,
                                                           0,
                                                           double,
                                                           2,
                                                           Dune::XT::LA::Backends::istl_dense>::bind(m);
  Dune::XT::Common::bindings::add_initialization(m, "dune.gdt.operators.elliptic");
}
