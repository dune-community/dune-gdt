# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
set(la_containerfactory_includedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/gdt/la/containerfactory)
set(la_containerfactory_include_HEADERS  common.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${la_containerfactory_include_HEADERS} DESTINATION ${la_containerfactory_includedir})