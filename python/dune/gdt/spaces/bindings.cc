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

#include <dune/pybindxi/pybind11.h>
#include <dune/pybindxi/stl.h>

#include <python/dune/xt/common/bindings.hh>
#include <dune/xt/grid/grids.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <python/dune/gdt/shared.hh>
#include <dune/gdt/spaces.hh>


#define DUNE_GDT_SPACES_BIND(_m, _GRID, _layer, _space_type, _space_backend, _p, _r, _rC, _grid_backend)               \
  Dune::GDT::bindings::SpaceInterface<Dune::GDT::SpaceProvider<_GRID,                                                  \
                                                               Dune::XT::Grid::Layers::_layer,                         \
                                                               Dune::GDT::SpaceType::_space_type,                      \
                                                               Dune::GDT::Backends::_space_backend,                    \
                                                               _p,                                                     \
                                                               double,                                                 \
                                                               _r,                                                     \
                                                               _rC,                                                    \
                                                               Dune::XT::Grid::Backends::_grid_backend>>::bind(_m)


PYBIND11_MODULE(__spaces, m)
{
  namespace py = pybind11;
  using namespace pybind11::literals;

  DUNE_GDT_SPACES_BIND(m, ALU_2D_SIMPLEX_CONFORMING, leaf, dg, gdt, 1, 1, 1, view);
  DUNE_GDT_SPACES_BIND(m, ALU_2D_SIMPLEX_CONFORMING, leaf, dg, gdt, 2, 1, 1, view);
  DUNE_GDT_SPACES_BIND(m, ALU_2D_SIMPLEX_CONFORMING, leaf, dg, gdt, 3, 1, 1, view);
  DUNE_GDT_SPACES_BIND(m, ALU_2D_SIMPLEX_CONFORMING, dd_subdomain, dg, gdt, 1, 1, 1, view);
  DUNE_GDT_SPACES_BIND(m, ALU_2D_SIMPLEX_CONFORMING, leaf, rt, gdt, 0, 2, 1, view);
  add_initialization(m, "dune.gdt.spaces");
}

#endif // HAVE_DUNE_PYBINDXI
