# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

enable_testing()

BEGIN_TESTCASES()

foreach (test ${testnames})
  add_dune_mpi_flags(${test})
  add_dune_alugrid_flags(${test})
  TARGET_LINK_LIBRARIES(${test} ${COMMON_LIBS})
endforeach (test ${testnames})

END_TESTCASES()

target_link_libraries(test_linearelliptic-cg-discretization_fem_eigen_alugrid    lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic-cg-discretization_fem_eigen_yaspgrid      lib_test_linearelliptic_cg_discretizations_yaspgrid)
target_link_libraries(test_linearelliptic-cg-discretization_fem_istl_alugrid     lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic-cg-discretization_fem_istl_yaspgrid       lib_test_linearelliptic_cg_discretizations_yaspgrid)
target_link_libraries(test_linearelliptic-cg-discretization_pdelab_eigen_alugrid lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic-cg-discretization_pdelab_eigen_yaspgrid   lib_test_linearelliptic_cg_discretizations_yaspgrid)
target_link_libraries(test_linearelliptic-cg-discretization_pdelab_istl_alugrid  lib_test_linearelliptic_cg_discretizations_alugrid)
target_link_libraries(test_linearelliptic-cg-discretization_pdelab_istl_yaspgrid    lib_test_linearelliptic_cg_discretizations_yaspgrid)

set_tests_properties(   test_linearelliptic-cg-discretization_pdelab_eigen_yaspgrid
                        test_linearelliptic-cg-discretization_pdelab_istl_yaspgrid
                        test_linearelliptic-cg-discretization_fem_eigen_yaspgrid
                        test_linearelliptic-cg-discretization_fem_istl_yaspgrid
                        test_linearelliptic-cg-discretization_fem_istl_alugrid
                        test_linearelliptic-cg-discretization_fem_eigen_alugrid
                        test_linearelliptic-cg-discretization_pdelab_eigen_alugrid
                        test_linearelliptic-cg-discretization_pdelab_istl_alugrid
                        test_elliptic-swipdg-discretization
    PROPERTIES TIMEOUT 720 )

# link spe10 data file if present
if (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
  set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
  if (NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
    execute_process(COMMAND ln -s "${SPE10MODEL1DATA}" "${SPE10MODEL1DATA_TARGET_FILENAME}")
  endif(NOT EXISTS ${SPE10MODEL1DATA_TARGET_FILENAME})
else()
    execute_process(COMMAND wget -nc http://wwwmath.uni-muenster.de/num/ohlberger/research/software/perm_case1.zip
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    execute_process(COMMAND unzip -o perm_case1.zip
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    set (SPE10MODEL1DATA_TARGET_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/perm_case1.dat")
endif (NOT ${SPE10MODEL1DATA} STREQUAL "SPE10MODEL1DATA-NOTFOUND")
