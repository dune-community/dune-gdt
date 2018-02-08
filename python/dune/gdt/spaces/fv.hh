// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef PYTHON_DUNE_GDT_SPACES_FV_BINDINGS_HH
#define PYTHON_DUNE_GDT_SPACES_FV_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <python/dune/xt/grid/grids.bindings.hh>

#include <dune/gdt/spaces/fv.hh>
#include <python/dune/gdt/spaces/interface.hh>


// begin: this is what we need for the .so

#define _DUNE_GDT_SPACES_FV_BIND_GDT(_m, _GRID, _layer, _r, _rC)                                                       \
  Dune::GDT::bindings::SpaceInterface<Dune::GDT::FvSpaceProvider<_GRID,                                                \
                                                                 Dune::XT::Grid::Layers::_layer,                       \
                                                                 Dune::GDT::Backends::gdt,                             \
                                                                 double,                                               \
                                                                 _r,                                                   \
                                                                 _rC>>::bind(_m)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA_LAYER(_m, _layer)                                                       \
//  _DUNE_GDT_SPACES_FV_BIND_GDT(_m, ALBERTA_2D, _layer, 1, 1)
//#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA(_m)                                                                     \
//  _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA_LAYER(_m, leaf);                                                              \
//  _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA_LAYER(_m, level)
//#else
#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA(_m)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALU_LAYER(_m, _layer)                                                             \
  _DUNE_GDT_SPACES_FV_BIND_GDT(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, 1, 1)
#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALU(_m)                                                                           \
  _DUNE_GDT_SPACES_FV_BIND_GDT_ALU_LAYER(_m, leaf);                                                                    \
  _DUNE_GDT_SPACES_FV_BIND_GDT_ALU_LAYER(_m, level)
#else
#define _DUNE_GDT_SPACES_FV_BIND_GDT_ALU(_m)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG
//#define _DUNE_GDT_SPACES_FV_BIND_GDT_UG_LAYER(_m, _layer) _DUNE_GDT_SPACES_FV_BIND_GDT(_m, UG_2D, _layer, 1, 1)
//#define _DUNE_GDT_SPACES_FV_BIND_GDT_UG(_m)                                                                          \
//  _DUNE_GDT_SPACES_FV_BIND_GDT_UG_LAYER(_m, leaf);                                                                   \
//  _DUNE_GDT_SPACES_FV_BIND_GDT_UG_LAYER(_m, level)
//#else
#define _DUNE_GDT_SPACES_FV_BIND_GDT_UG(_m)
//#endif

#define _DUNE_GDT_SPACES_FV_BIND_GDT_YASP_LAYER(_m, _layer)                                                            \
  _DUNE_GDT_SPACES_FV_BIND_GDT(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, 1, 1)
#define _DUNE_GDT_SPACES_FV_BIND_GDT_YASP(_m)                                                                          \
  _DUNE_GDT_SPACES_FV_BIND_GDT_YASP_LAYER(_m, leaf);                                                                   \
  _DUNE_GDT_SPACES_FV_BIND_GDT_YASP_LAYER(_m, level)

#define DUNE_GDT_SPACES_FV_BIND(_m)                                                                                    \
  _DUNE_GDT_SPACES_FV_BIND_GDT_ALBERTA(_m);                                                                            \
  _DUNE_GDT_SPACES_FV_BIND_GDT_ALU(_m);                                                                                \
  _DUNE_GDT_SPACES_FV_BIND_GDT_UG(_m);                                                                                 \
  _DUNE_GDT_SPACES_FV_BIND_GDT_YASP(_m)

// end: this is what we need for the .so


#endif // HAVE_DUNE_PYBINDXI
#endif // PYTHON_DUNE_GDT_SPACES_FV_BINDINGS_HH
