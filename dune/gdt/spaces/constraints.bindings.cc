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

#include "constraints.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .cc source file
DUNE_GDT_SPACES_CONSTRAINTS_BIND(template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view);

DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, common_dense);
#if HAVE_EIGEN
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, eigen_dense);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, eigen_sparse);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, YASP_2D_EQUIDISTANT_OFFSET, leaf, view, istl_sparse);
#endif

#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_FEM(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET);
#endif


#if HAVE_DUNE_ALUGRID
DUNE_GDT_SPACES_CONSTRAINTS_BIND(template, ALU_2D_SIMPLEX_CONFORMING, leaf, view);
DUNE_GDT_SPACES_CONSTRAINTS_BIND(template, ALU_2D_SIMPLEX_CONFORMING, level, view);

DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, common_dense);
#if HAVE_EIGEN
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, eigen_dense);
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, eigen_sparse);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_LA(template, ALU_2D_SIMPLEX_CONFORMING, leaf, view, istl_sparse);
#endif

#if HAVE_DUNE_FEM
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_FEM(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#if HAVE_DUNE_PDELAB
DUNE_GDT_SPACES_CONSTRAINTS_ADDBIND_ASSEMBLER_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING);
#endif
#endif


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
