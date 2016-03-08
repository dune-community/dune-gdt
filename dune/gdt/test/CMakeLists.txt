# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

enable_testing()

BEGIN_TESTCASES()

END_TESTCASES()

foreach (test ${testnames})
  add_dune_alugrid_flags(${test})
  TARGET_LINK_LIBRARIES(${test} ${COMMON_LIBS})
  set_tests_properties(${test} PROPERTIES TIMEOUT 360)
  if (${test} MATCHES spaces)
    list(APPEND spaces_test_binaries ${test})
  endif (${test})
  if (${test} MATCHES operators)
    list(APPEND operators_test_binaries ${test})
  endif (${test})
  if (${test} MATCHES projections)
    list(APPEND projections_test_binaries ${test})
  endif (${test})
  if (${test} MATCHES prolongations)
    list(APPEND prolongations_test_binaries ${test})
  endif (${test})
  if (${test} MATCHES linearelliptic AND ${test} MATCHES alugrid)
    list(APPEND linearelliptic_alugrid_test_binaries ${test})
  endif (${test})
  if (${test} MATCHES linearelliptic AND ${test} MATCHES sgrid)
    list(APPEND linearelliptic_sgrid_test_binaries ${test})
  endif (${test})
endforeach (test ${testnames})

add_custom_target(test_binaries_space DEPENDS ${spaces_test_binaries})
add_custom_target(test_binaries_operator DEPENDS ${operators_test_binaries})
add_custom_target(test_binaries_projection DEPENDS ${projections_test_binaries})
add_custom_target(test_binaries_prolongation DEPENDS ${prolongations_test_binaries})
add_custom_target(test_binaries_linearelliptic_cg_discretization_alugrid DEPENDS ${linearelliptic_alugrid_test_binaries})
add_custom_target(test_binaries_linearelliptic_cg_discretization_sgrid DEPENDS ${linearelliptic_sgrid_test_binaries})

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")

add_subdirectory(linearelliptic)
add_subdirectory(hyperbolic)

target_link_libraries(test_linearelliptic__cg_discretization__alugrid_fem_eigen    lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic__cg_discretization__sgrid_fem_eigen      lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic__cg_discretization__alugrid_fem_istl     lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic__cg_discretization__sgrid_fem_istl       lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic__cg_discretization__alugrid_pdelab_eigen lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic__cg_discretization__sgrid_pdelab_eigen   lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic__cg_discretization__alugrid_pdelab_istl  lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic__cg_discretization__sgrid_pdelab_istl    lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_hyperbolic__fv_discretization__yaspgrid                 lib_test_hyperbolic_fv_discretizations_yaspgrid)
