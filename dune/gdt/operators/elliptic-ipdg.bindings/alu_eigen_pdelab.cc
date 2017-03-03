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
#if HAVE_DUNE_ALUGRID
#if HAVE_EIGEN
#if HAVE_DUNE_PDELAB

#include "../elliptic-ipdg.bindings.hh"

namespace Dune {
namespace GDT {
namespace bindings {


DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING, eigen_dense);
DUNE_GDT_OPERATORS_ELLIPTIC_IPDG_BIND_PDELAB(template, ALU_2D_SIMPLEX_CONFORMING, eigen_sparse);


} // namespace bindings
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_PDELAB
#endif // HAVE_EIGEN
#endif // HAVE_DUNE_ALUGRID
#endif // HAVE_DUNE_PYBINDXI
