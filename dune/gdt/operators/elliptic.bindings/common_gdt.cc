// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include "../elliptic.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET, common_dense);

#if HAVE_DUNE_ALUGRID
DUNE_GDT_OPERATORS_ELLIPTIC_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING, common_dense);
#endif // HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
