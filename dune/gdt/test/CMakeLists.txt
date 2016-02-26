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
endforeach (test ${testnames})

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")

add_subdirectory(linearelliptic)

target_link_libraries(test_linearelliptic_cg_discretization_fem_eigen_alugrid    lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic_cg_discretization_fem_eigen_sgrid      lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic_cg_discretization_fem_istl_alugrid     lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic_cg_discretization_fem_istl_sgrid       lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic_cg_discretization_pdelab_eigen_alugrid lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic_cg_discretization_pdelab_eigen_sgrid   lib_test_linearelliptic_cg_discretizations_sgrid)
target_link_libraries(test_linearelliptic_cg_discretization_pdelab_istl_alugrid  lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic_cg_discretization_pdelab_istl_sgrid    lib_test_linearelliptic_cg_discretizations_sgrid)

