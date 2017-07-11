// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#include <config.h>

#if HAVE_LAPACK

#include "eigensolver.hh"

#include <lapacke.h>

namespace Dune {
namespace GDT {


int LapackWrapper::dggev(char jobvl,
                         char jobvr,
                         int n,
                         double* a,
                         int lda,
                         double* b,
                         int ldb,
                         double* alphar,
                         double* alphai,
                         double* beta,
                         double* vl,
                         int ldvl,
                         double* vr,
                         int ldvr)
{
  return LAPACKE_dggev(LAPACK_ROW_MAJOR, jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr);
}


} // namespace GDT
} // namespace Dune

#endif // HAVE_LAPACK
