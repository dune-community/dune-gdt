# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2014, 2016)

set(M4FILES  dune-gdt.m4)
set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/aclocal)
set(aclocal_DATA  ${M4FILES})
set(EXTRA_DIST  ${M4FILES})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})