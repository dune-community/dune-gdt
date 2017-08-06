# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013, 2016 - 2017)
#   Rene Milk       (2016 - 2017)
#   Tobias Leibner  (2016)

enable_testing()

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  BEGIN_TESTCASES(dunegdt)
else()
  BEGIN_TESTCASES()
endif()

END_TESTCASES()

# (re)generated from refresh_test_timings target
if(DEFINED ENV{TRAVIS})
    include("builder_definitions.cmake")
    # make sure there's the same amount of builders in travis as bins from distribution
    add_custom_target(regenerate_travis_yml ${CMAKE_BINARY_DIR}/dune-env
                        "${CMAKE_SOURCE_DIR}/.travis.yml.py" "${DXT_BIN_COUNT}" )
    add_dependencies(refresh_test_timings regenerate_travis_yml)
endif(DEFINED ENV{TRAVIS})

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
    if(${test_target} MATCHES linearelliptic__block_swipdg_discretization)
        target_link_libraries(${test_target} block_swipdg_expectations)
    endif()
    if(${test_target} MATCHES linearelliptic__cg_discretization)
        target_link_libraries(${test_target} cg_expectations)
    endif()
    if(${test_target} MATCHES linearelliptic__swipdg_discretization)
        target_link_libraries(${test_target} swipdg_expectations)
    endif()
    if(${test_target} MATCHES test_linearelliptic__swipdg_estimators)
        target_link_libraries(${test_target} swipdg_expectations)
        target_link_libraries(${test_target} swipdg_estimator_expectations)
    endif()
    if(${test_target} MATCHES hyperbolic__fv_discretization)
	    target_link_libraries(${test_target} ${LPSolve_LIBRARIES} ${Qhull_LIBRARIES} ${LAPACKE_LIBRARIES} ${MATEXP_LIBRARIES} ${LEBEDEV_DATA_LIBRARIES} expectations_hyperbolic_fv_discretizations_yaspgrid)
    endif()
    if(${test_target} MATCHES entropytest)
	    target_link_libraries(${test_target} ${LPSolve_LIBRARIES} ${Qhull_LIBRARIES} ${LAPACKE_LIBRARIES} ${MATEXP_LIBRARIES} ${LEBEDEV_DATA_LIBRARIES})
    endif()
endforeach()
