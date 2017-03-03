// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_FV_BINDINGS_HH
#define DUNE_GDT_SPACES_FV_BINDINGS_HH
#if HAVE_DUNE_PYBINDXI

#include <dune/xt/grid/grids.bindings.hh>

#include "fv.hh"
#include "interface.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


#define _DEFINE_FV_SPACE(_G, _layer, _backend, _r, _rC)                                                                \
  typedef typename FvSpaceProvider<_G, XT::Grid::Layers::_layer, ChooseSpaceBackend::_backend, double, _r, _rC>::type  \
      Fv_##_G##_##_layer##_to_##_r##x##_rC##_##_backend##_Space

_DEFINE_FV_SPACE(YASP_2D_EQUIDISTANT_OFFSET, leaf, gdt, 1, 1);
_DEFINE_FV_SPACE(YASP_2D_EQUIDISTANT_OFFSET, leaf, gdt, 1, 1);
_DEFINE_FV_SPACE(YASP_2D_EQUIDISTANT_OFFSET, level, gdt, 1, 1);
_DEFINE_FV_SPACE(YASP_2D_EQUIDISTANT_OFFSET, level, gdt, 1, 1);
#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
_DEFINE_FV_SPACE(ALU_2D_SIMPLEX_CONFORMING, leaf, gdt, 1, 1);
_DEFINE_FV_SPACE(ALU_2D_SIMPLEX_CONFORMING, leaf, gdt, 1, 1);
_DEFINE_FV_SPACE(ALU_2D_SIMPLEX_CONFORMING, level, gdt, 1, 1);
_DEFINE_FV_SPACE(ALU_2D_SIMPLEX_CONFORMING, level, gdt, 1, 1);
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID

#undef _DEFINE_FV_SPACE


// this is used by other headers
#define FV_SPACE(_G, _layer, _backend, _r, _rC) Fv_##_G##_##_layer##_to_##_r##x##_rC##_##_backend##_Space


#define DUNE_GDT_SPACES_FV_BIND_GDT(_prefix, _GRID)                                                                    \
  _prefix class SpaceInterface<FV_SPACE(_GRID, leaf, gdt, 1, 1)>;                                                      \
  _prefix class SpaceInterface<FV_SPACE(_GRID, level, gdt, 1, 1)>


// these lines have to match the corresponding ones in the .cc source file
DUNE_GDT_SPACES_FV_BIND_GDT(extern template, YASP_2D_EQUIDISTANT_OFFSET);

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_SPACES_FV_BIND_GDT(extern template, ALU_2D_SIMPLEX_CONFORMING);
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
#endif // DUNE_GDT_SPACES_FV_BINDINGS_HH
