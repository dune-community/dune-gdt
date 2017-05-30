// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_RT_BINDINGS_HH
#define DUNE_GDT_SPACES_RT_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/xt/grid/grids.bindings.hh>

#include "interface.bindings.hh"
#include "rt.hh"


// begin: this is what we need for the .so

#if HAVE_DUNE_PDELAB
#define _DUNE_GDT_SPACES_RT_BIND_PDELAB(_m, _GRID, _layer)                                                             \
  Dune::GDT::bindings::SpaceInterface<Dune::GDT::RtSpaceProvider<_GRID,                                                \
                                                                 Dune::XT::Grid::Layers::_layer,                       \
                                                                 Dune::GDT::Backends::pdelab,                          \
                                                                 0,                                                    \
                                                                 double,                                               \
                                                                 _GRID::dimension,                                     \
                                                                 1>>::bind(_m)

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU_LAYER(_m, _layer)                                                          \
  _DUNE_GDT_SPACES_RT_BIND_PDELAB(_m, ALU_2D_SIMPLEX_CONFORMING, _layer)
#define _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU(_m)                                                                        \
  _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU_LAYER(_m, leaf);                                                                 \
  _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU_LAYER(_m, level)
#else
#define _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU(_m)
#endif

#define _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALL(_m) _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALU(_m)
#else // HAVE_DUNE_PDELAB
#define _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALL(_m)
#endif

#define DUNE_GDT_SPACES_RT_BIND(_m) _DUNE_GDT_SPACES_RT_BIND_PDELAB_ALL(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_RT_BINDINGS_HH
