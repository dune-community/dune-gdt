// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#include <dune/gdt/functionals/elliptic-ipdg.bindings.hh>


#if HAVE_DUNE_ALUGRID && HAVE_DUNE_FEM && HAVE_DUNE_ISTL
DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(template, leaf, part, dg, fem, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(template, level, part, dg, fem, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_ELLIPTIC_IPDG_BIND_LIB_ALU(template, dd_subdomain, part, dg, fem, 1, istl_sparse);
#endif


#endif // HAVE_DUNE_PYBINDXI
