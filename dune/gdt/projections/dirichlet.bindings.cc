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

#include "dirichlet.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


// these lines have to match the corresponding ones in the .hh header file
#if HAVE_DUNE_FEM
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, YASP_2D_EQUIDISTANT_OFFSET, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB

#if HAVE_DUNE_ALUGRID
#if HAVE_DUNE_FEM
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_FEM(template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_FEM
#if HAVE_DUNE_PDELAB
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING, COMMON_DENSE_VECTOR);
#if HAVE_EIGEN
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING, EIGEN_DENSE_VECTOR);
#endif
#if HAVE_DUNE_ISTL
DUNE_GDT_PROJECTIONS_DIRICHLET_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING, ISTL_DENSE_VECTOR);
#endif
#endif // HAVE_DUNE_PDELAB
#endif // HAVE_DUNE_ALUGRID


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PYBINDXI
