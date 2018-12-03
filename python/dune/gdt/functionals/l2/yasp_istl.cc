// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017 - 2018)
//   René Fritze     (2018)

#include "config.h"

#if HAVE_DUNE_PYBINDXI

#  include <python/dune/gdt/functionals/l2/bindings.hh>


#  if HAVE_DUNE_ISTL
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(template, leaf, view, cg, gdt, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(template, level, view, cg, gdt, 1, istl_sparse);
DUNE_GDT_FUNCTIONALS_L2_BIND_LIB_YASP(template, dd_subdomain, view, cg, gdt, 1, istl_sparse);
#  endif


#endif // HAVE_DUNE_PYBINDXI
