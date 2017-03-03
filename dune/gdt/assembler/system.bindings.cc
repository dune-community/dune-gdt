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

#include "system.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(template, YASP_2D_EQUIDISTANT_OFFSET);
#if HAVE_DUNE_FEM
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_GDT(template, ALU_2D_SIMPLEX_CONFORMING);
#if HAVE_DUNE_FEM
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_ASSEMBLER_SYSTEM_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif // HAVE_ALUGRID || HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
