# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2013, 2016 - 2017)
#   Rene Milk       (2016 - 2018)
#   Tobias Leibner  (2016 - 2017)
# ~~~

enable_testing()

set(DXT_TEST_TIMEOUT 1000)
if(DUNE_XT_WITH_PYTHON_BINDINGS)
  begin_testcases(dunegdt)
else()
  begin_testcases()
endif()

end_testcases()

# link spe10 data file if present
if(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set(SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
else(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  message(STATUS "Downloading Spe10 model1 data file:" ${SPE10MODEL1DATA})
  execute_process(COMMAND wget -nc https://www.uni-muenster.de/AMM/mirror/spe10/model1/perm_case1.zip
                  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  execute_process(COMMAND unzip -o perm_case1.zip WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set(SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
endif(NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")

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
  if(${test_target} MATCHES linearelliptic__swipdg_estimators)
    target_link_libraries(${test_target} swipdg_expectations)
    target_link_libraries(${test_target} swipdg_estimator_expectations)
  endif()
  if(${test_target} MATCHES hyperbolic__fv_discretization)
    target_link_libraries(${test_target}
                          ${Qhull_LIBRARIES}
                          ${LAPACKE_LIBRARIES}
                          ${MATEXP_LIBRARIES}
                          ${LEBEDEVDATA_LIBRARIES}
                          expectations_hyperbolic_fv_discretizations_yaspgrid)
  endif()
  if(${test_target} MATCHES hyperbolic__momentmodels)
    target_link_libraries(${test_target}
                          ${Qhull_LIBRARIES}
                          ${LAPACKE_LIBRARIES}
                          ${MATEXP_LIBRARIES}
                          ${FEKETE_LIBRARIES}
                          ${LEBEDEVDATA_LIBRARIES})
  endif()
  if(${test_target} MATCHES entropytest)
    target_link_libraries(${test_target}
                          ${Qhull_LIBRARIES}
                          ${LAPACKE_LIBRARIES}
                          ${MATEXP_LIBRARIES}
                          ${FEKETE_LIBRARIES}
                          ${LEBEDEVDATA_LIBRARIES})
  endif()
endforeach()
