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

#include <dune/xt/common/bindings.hh>

#include "block.bindings.hh"


#define DUNE_GDT_SPACES_BLOCK_BIND(_m, _GRID, _s_type, _s_backend, _p)                                                 \
  Dune::GDT::bindings::BlockSpace<Dune::GDT::SpaceProvider<_GRID,                                                      \
                                                           Dune::XT::Grid::Layers::dd_subdomain,                       \
                                                           Dune::GDT::SpaceType::_s_type,                              \
                                                           Dune::GDT::Backends::_s_backend,                            \
                                                           _p,                                                         \
                                                           double,                                                     \
                                                           1,                                                          \
                                                           1>>::bind(_m)


PYBIND11_PLUGIN(__spaces_block)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  py::module m("__spaces_block", "dune-gdt: Block spaces");
  DUNE_XT_COMMON_BINDINGS_INITIALIZE(m, "dune.gdt.spaces.block");

  DUNE_GDT_SPACES_BLOCK_BIND(m, ALU_2D_SIMPLEX_CONFORMING, dg, fem, 1);

  return m.ptr();
}

#endif // HAVE_DUNE_PYBINDXI
