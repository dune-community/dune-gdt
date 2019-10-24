# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013, 2016 - 2018)
#   René Fritze     (2016 - 2018)
#   Tobias Leibner  (2016 - 2018)
# ~~~

enable_testing()

set(DXT_TEST_TIMEOUT 1000)

add_subdir_tests(burgers)
add_subdir_tests(discretefunction)
add_subdir_tests(instationary-eocstudies)
add_subdir_tests(integrands)
add_subdir_tests(interpolations)
add_subdir_tests(inviscid-compressible-flow)
add_subdir_tests(linear-transport)
add_subdir_tests(misc)
add_subdir_tests(momentmodels)
add_subdir_tests(operators)
add_subdir_tests(projections)
add_subdir_tests(prolongations)
add_subdir_tests(spaces)
add_subdir_tests(stationary-eocstudies)
add_subdir_tests(stationary-heat-equation)
add_subdir_tests(stokes)

finalize_test_setup()
