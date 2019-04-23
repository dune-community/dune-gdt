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

add_library("boltzmann" SHARED ${lib_boltzmann_sources})
target_link_libraries(boltzmann ${DUNE_LIBS})
add_dune_mpi_flags(boltzmann)
if(TBB_FOUND)
  add_dune_tbb_flags(boltzmann)
endif(TBB_FOUND)

if(dune-grid_FOUND)
  add_dune_superlu_flags(boltzmann)
endif(dune-grid_FOUND)
target_link_dune_default_libraries(boltzmann)

find_package(PythonLibs 3)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(boltzmann
                      ${PYTHON_LIBRARIES}
                      ${Qhull_LIBRARIES}
                      ${LAPACKE_LIBRARIES}
                      ${MATEXP_LIBRARIES}
                      ${FEKETE_LIBRARIES}
                      ${LEBEDEVDATA_LIBRARIES})
