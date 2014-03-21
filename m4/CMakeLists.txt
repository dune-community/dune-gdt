# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
set(M4FILES  dune-gdt.m4)
set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/aclocal)
set(aclocal_DATA  ${M4FILES})
set(EXTRA_DIST  ${M4FILES})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})