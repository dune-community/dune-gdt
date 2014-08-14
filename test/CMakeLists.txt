# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

enable_testing()

#HEADERCHECK( ${gdt} )
#DEPENDENCYCHECK( ${gdt} )

BEGIN_TESTCASES()

foreach( test ${testnames} )
  IF(NOT CMAKE_WITH_AUTOTOOLS)
    add_dune_mpi_flags(${test})
    add_dune_alugrid_flags(${test})
  ENDIF(NOT CMAKE_WITH_AUTOTOOLS)
  TARGET_LINK_LIBRARIES(${test} ${COMMON_LIBS})
endforeach( test ${testnames} )

END_TESTCASES()
