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

#include <dune/common/parallel/mpihelper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/misc/mpimanager.hh>
#endif

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <dune/xt/grid/dd/subdomains/grid.hh>
#include <dune/xt/grid/grids.hh>

#include <dune/gdt/spaces.hh>
#include <python/dune/gdt/operators/fluxreconstruction.bindings.hh>


PYBIND11_PLUGIN(__operators_fluxreconstruction)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__operators_fluxreconstruction", "dune-gdt: DiffusiveFluxReconstructionOperator");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.operators.fluxreconstruction");

#if HAVE_DUNE_ALUGRID && HAVE_DUNE_ISTL
  Dune::GDT::bindings::DiffusiveFluxReconstructionOperator<ALU_2D_SIMPLEX_CONFORMING,
                                                           Dune::GDT::SpaceType::rt,
                                                           Dune::GDT::Backends::gdt,
                                                           Dune::XT::Grid::Layers::leaf,
                                                           0,
                                                           double,
                                                           2,
                                                           Dune::XT::LA::Backends::istl_dense>::bind(m);
#endif // HAVE_DUNE_ALUGRID && HAVE_DUNE_ISTL

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
