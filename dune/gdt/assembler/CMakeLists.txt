# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
set(assembler_includedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/gdt/assembler)
set(assembler_include_HEADERS  system.hh)
set(SUBDIRS  local)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${assembler_include_HEADERS} DESTINATION ${assembler_includedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
  unset(opt)
endforeach(i ${SUBDIRS})
