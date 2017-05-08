// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI && HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM

#include "../system.bindings.hh"


namespace Dune {
namespace GDT {
namespace bindings {


DUNE_GDT_ASSEMBLER_SYSTEM_BIND_LIB_ALU_FEM(template);


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI && HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM
