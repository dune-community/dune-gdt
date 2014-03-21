# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
set(discretefunction_includedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/gdt/discretefunction)
set(discretefunction_include_HEADERS  default.hh local.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${discretefunction_include_HEADERS} DESTINATION ${discretefunction_includedir})