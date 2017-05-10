// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB_HH
#define DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB_HH

#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/layers.hh>

#include "dune-fem-wrapper.hh"

#if DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_DUNE_FEM


#define _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, _layer_type, _p, _R, _r, _rC)                         \
  _prefix class Dune::GDT::DuneFemCgSpaceWrapper<                                                                      \
      typename Dune::XT::Grid::Layer<_GRID,                                                                            \
                                     Dune::XT::Grid::Layers::_layer_type,                                              \
                                     Dune::XT::Grid::Backends::part,                                                   \
                                     Dune::XT::Grid::DD::SubdomainGrid<_GRID>>::type,                                  \
      _p,                                                                                                              \
      _R,                                                                                                              \
      _r,                                                                                                              \
      _rC>

#define _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB_ALL_LAYERS(_prefix, _GRID, _p, _R, _r, _rC)                           \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, adaptive_leaf, _p, _R, _r, _rC);                            \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, leaf, _p, _R, _r, _rC);                                     \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, level, _p, _R, _r, _rC);                                    \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, dd_subdomain, _p, _R, _r, _rC);                             \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, dd_subdomain_boundary, _p, _R, _r, _rC);                    \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID, dd_subdomain_coupling, _p, _R, _r, _rC)

#define DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(_prefix, _GRID)                                                        \
  _DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB_ALL_LAYERS(_prefix, _GRID, 1, double, 1, 1)

#if HAVE_DUNE_ALUGRID
DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif
DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(extern template, YASP_1D_EQUIDISTANT_OFFSET);
DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(extern template, YASP_2D_EQUIDISTANT_OFFSET);
DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB(extern template, YASP_3D_EQUIDISTANT_OFFSET);


#endif // DUNE_XT_WITH_PYTHON_BINDINGS && HAVE_DUNE_FEM


#endif // DUNE_GDT_SPACES_CG_DUNE_FEM_WRAPPER_LIB_HH
