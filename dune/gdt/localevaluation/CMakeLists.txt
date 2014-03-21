# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
set(includedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/gdt/evaluation/local)
set(include_HEADERS  elliptic.hh interface.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${include_HEADERS} DESTINATION ${includedir})