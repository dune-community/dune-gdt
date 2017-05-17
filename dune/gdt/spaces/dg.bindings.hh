// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_DG_BINDINGS_HH
#define DUNE_GDT_SPACES_DG_BINDINGS_HH
//#if HAVE_DUNE_PYBINDXI

#include <dune/xt/grid/grids.bindings.hh>

#include "dg.hh"
#include "interface.bindings.hh"


// begin: this is what we need for the .so

// * fem
#if HAVE_DUNE_FEM
#define _DUNE_GDT_SPACES_DG_BIND_FEM(_m, _GRID, _layer, _r, _rC)                                                       \
  Dune::GDT::bindings::SpaceInterface<Dune::GDT::DgSpaceProvider<_GRID,                                                \
                                                                 Dune::XT::Grid::Layers::_layer,                       \
                                                                 Dune::GDT::Backends::fem,                             \
                                                                 1,                                                    \
                                                                 double,                                               \
                                                                 _r,                                                   \
                                                                 _rC>>::bind(_m)

//#if HAVE_ALBERTA
//#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA_LAYER(_m, _layer)                                                       \
//  _DUNE_GDT_SPACES_DG_BIND_FEM(_m, ALBERTA_2D, _layer, 1, 1)
//#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA(_m)                                                                     \
//  _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA_LAYER(_m, dd_subdomain);                                                      \
//  _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA_LAYER(_m, leaf);                                                              \
//  _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA_LAYER(_m, level)
//#else
#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA(_m)
//#endif

#if HAVE_DUNE_ALUGRID
#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALU_LAYER(_m, _layer)                                                             \
  _DUNE_GDT_SPACES_DG_BIND_FEM(_m, ALU_2D_SIMPLEX_CONFORMING, _layer, 1, 1)
#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALU(_m)                                                                           \
  _DUNE_GDT_SPACES_DG_BIND_FEM_ALU_LAYER(_m, dd_subdomain);                                                            \
  _DUNE_GDT_SPACES_DG_BIND_FEM_ALU_LAYER(_m, leaf);                                                                    \
  _DUNE_GDT_SPACES_DG_BIND_FEM_ALU_LAYER(_m, level)
#else
#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALU(_m)
#endif

//#if HAVE_DUNE_UGGRID || HAVE_UG // <- does not work
//#define _DUNE_GDT_SPACES_DG_BIND_FEM_UG_LAYER(_m, _layer) _DUNE_GDT_SPACES_DG_BIND_FEM(_m, UG_2D, _layer, 1, 1)
//#define _DUNE_GDT_SPACES_DG_BIND_FEM_UG(_m)
//  _DUNE_GDT_SPACES_DG_BIND_FEM_UG_LAYER(_m, dd_subdomain);
//  _DUNE_GDT_SPACES_DG_BIND_FEM_UG_LAYER(_m, leaf);
//  _DUNE_GDT_SPACES_DG_BIND_FEM_UG_LAYER(_m, level)
//#else
//#define _DUNE_GDT_SPACES_DG_BIND_FEM_UG(_m)
//#endif

#define _DUNE_GDT_SPACES_DG_BIND_FEM_YASP_LAYER(_m, _layer)                                                            \
  _DUNE_GDT_SPACES_DG_BIND_FEM(_m, YASP_1D_EQUIDISTANT_OFFSET, _layer, 1, 1);                                          \
  _DUNE_GDT_SPACES_DG_BIND_FEM(_m, YASP_2D_EQUIDISTANT_OFFSET, _layer, 1, 1)
#define _DUNE_GDT_SPACES_DG_BIND_FEM_YASP(_m)
//  _DUNE_GDT_SPACES_DG_BIND_FEM_YASP_LAYER(_m, dd_subdomain);                                                           \
//  _DUNE_GDT_SPACES_DG_BIND_FEM_YASP_LAYER(_m, leaf);                                                                   \
//  _DUNE_GDT_SPACES_DG_BIND_FEM_YASP_LAYER(_m, level)

#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALL(_m)                                                                           \
  _DUNE_GDT_SPACES_DG_BIND_FEM_ALBERTA(_m);                                                                            \
  _DUNE_GDT_SPACES_DG_BIND_FEM_ALU(_m);                                                                                \
  _DUNE_GDT_SPACES_DG_BIND_FEM_YASP(_m)
//_DUNE_GDT_SPACES_DG_BIND_FEM_UG(_m); // <- does not work
#else // HAVE_DUNE_FEM
#define _DUNE_GDT_SPACES_DG_BIND_FEM_ALL(_m)
#endif

#define DUNE_GDT_SPACES_DG_BIND(_m) _DUNE_GDT_SPACES_DG_BIND_FEM_ALL(_m)

// end: this is what we need for the .so


//#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_DG_BINDINGS_HH
