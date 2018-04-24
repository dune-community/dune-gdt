# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Rene Milk (2017 - 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

Comparison of testcase properties
=================================


|                     | AO2013 | ER2007 | ESV2007 | Spe10Model1 | MixedBoundary |
|---------------------|:------:|:------:|:-------:|:-----------:|:-------------:|
| UnitCube            |   ✓    |   ✓    |    x    |      x      |       ✓       |
| Constant Diffusion  |   x    |   ✓    |    ✓    |      x      |       ✓       |
| non-zero dirichlet  |   x    |   ✓    |    x    |      x      |       ✓       |
| dirichlet only      |   ✓    |   ✓    |    ✓    |      ✓      |       x       |
| exact solution      |   x    |   ✓    |    ✓    |      x      |       x       |
| DD testcase         |   x    |   ✓    |    ✓    |      x      |       x       |
