# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2013, 2016)
#   Tobias Leibner  (2016)

enable_testing()

BEGIN_TESTCASES()

END_TESTCASES()

if ($ENV{TRAVIS})
  # (re)generated from refresh_test_timings target
  include("builder_definitions.cmake")
endif ()

# make sure there's the same amount of builders in travis as bins from distribution
add_custom_target(regenerate_travis_yml ${CMAKE_BINARY_DIR}/dune-env-3
                    "${CMAKE_SOURCE_DIR}/.travis.yml.py" "${DXT_BIN_COUNT}" )
add_dependencies(refresh_test_timings regenerate_travis_yml)

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
else(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
        message(STATUS "Downloading Spe10 model1 data file:" ${SPE10MODEL1DATA})
  execute_process(COMMAND wget -nc http://wwwmath.uni-muenster.de/num/ohlberger/research/software/perm_case1.zip
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND unzip -o perm_case1.zip WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")

add_subdirectory(linearelliptic)
add_subdirectory(hyperbolic)


foreach(test_target ${dxt_test_binaries})
    if(${test_target} MATCHES linearelliptic__)
        target_link_libraries(${test_target} cg_expectations)
    endif(${test_target} MATCHES linearelliptic__)
    if(${test_target} MATCHES linearelliptic__)
        target_link_libraries(${test_target} swipdg_expectations)
    endif(${test_target} MATCHES linearelliptic__)
    if(${test_target} MATCHES hyperbolic__fv_discretization)
        target_link_libraries(${test_target} expectations_hyperbolic_fv_discretizations_yaspgrid)
    endif(${test_target} MATCHES hyperbolic__fv_discretization)
endforeach(test_target ${dxt_test_binaries})

