# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2010, 2013 - 2017)
#   Rene Milk       (2017 - 2018)
#   Tobias Leibner  (2017)

if(DUNE_XT_WITH_PYTHON_BINDINGS)
  set(dunegdt_sources
      assembler/system.lib/alu_2d_simplex_conforming.cc
      assembler/system.lib/yasp_1d_equidistant_offset.cc
      assembler/system.lib/yasp_2d_equidistant_offset.cc
      assembler/system.lib/yasp_3d_equidistant_offset.cc
      spaces/cg.lib/alu_2d_simplex_conforming.cc
      spaces/cg.lib/yasp_1d_equidistant_offset.cc
      spaces/cg.lib/yasp_2d_equidistant_offset.cc
      spaces/cg.lib/yasp_3d_equidistant_offset.cc
      spaces/dg.lib/alu_2d_simplex_conforming.cc
      spaces/dg.lib/yasp_1d_equidistant_offset.cc
      spaces/dg.lib/yasp_2d_equidistant_offset.cc
      spaces/dg.lib/yasp_3d_equidistant_offset.cc)
  dune_library_add_sources(dunegdt SOURCES ${dunegdt_sources})
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)

