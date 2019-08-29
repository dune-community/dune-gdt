# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2018)
#   René Fritze     (2017 - 2018)
#   Tobias Leibner  (2017 - 2018)
# ~~~

add_subdirectory(test EXCLUDE_FROM_ALL)

set(lib_boltzmann_sources test/boltzmanntest.cc)

dune_library_add_sources("boltzmann" SOURCES ${lib_boltzmann_sources})
