// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk (2018)

#ifndef PYTHON_DUNE_GDT_SPACES_RT_BINDINGS_HH
#define PYTHON_DUNE_GDT_SPACES_RT_BINDINGS_HH

#include <python/dune/xt/grid/grids.bindings.hh>

#include <python/dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/rt.hh>


// begin: this is what we need for the .so

#define _DUNE_GDT_SPACES_RT_BIND(_m, _GRID, _layer)                                                                    \
  Dune::GDT::bindings::SpaceInterface<Dune::GDT::RtSpaceProvider<_GRID,                                                \
                                                                 Dune::XT::Grid::Layers::_layer,                       \
                                                                 Dune::GDT::Backends::gdt,                             \
                                                                 0,                                                    \
                                                                 double,                                               \
                                                                 _GRID::dimension,                                     \
                                                                 1>>::bind(_m)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_SPACES_RT_BIND_ALU_LAYER(_m, _layer) _DUNE_GDT_SPACES_RT_BIND(_m, ALU_2D_SIMPLEX_CONFORMING, _layer); \
_DUNE_GDT_SPACES_RT_BIND(_m, ALU_2D_SIMPLEX_NONCONFORMING, _layer)
#define _DUNE_GDT_SPACES_RT_BIND_ALU(_m)                                                                               \
  _DUNE_GDT_SPACES_RT_BIND_ALU_LAYER(_m, leaf);                                                                        \
  _DUNE_GDT_SPACES_RT_BIND_ALU_LAYER(_m, level)
#else
#define _DUNE_GDT_SPACES_RT_BIND_ALU(_m)
#endif

#define DUNE_GDT_SPACES_RT_BIND(_m) _DUNE_GDT_SPACES_RT_BIND_ALU(_m)

// end: this is what we need for the .so


#endif // PYTHON_DUNE_GDT_SPACES_RT_BINDINGS_HH
