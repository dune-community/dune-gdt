// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2019)

#include "config.h"

#include <dune/common/parallel/mpihelper.hh>

#include <examples/coordinate-transformed-mn.hh>

int main(int argc, char* argv[])
{
  CoordinateTransformedBoltzmannSolver<3> solver;
  solver.solve();
} // ... main(...)
