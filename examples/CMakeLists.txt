# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2013, 2018)
#   René Fritze     (2018)
#   Tobias Leibner  (2018)
# ~~~

# Removes directory path, .cc suffix and leading __, replaces / with __
macro(CLEANUP_TARGET_NAME arg)
  string(REPLACE ${PROJECT_SOURCE_DIR} "" _targname ${arg})
  string(REPLACE ".cc" "" _targname ${_targname})
  string(REGEX REPLACE "/" "__" _targname ${_targname})
  string(REGEX REPLACE "^__" "" _targname ${_targname})
endmacro(CLEANUP_TARGET_NAME)

file(GLOB_RECURSE example_sources "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
foreach(source ${example_sources}) # do some name conversion
  set(_targname ${source})
  cleanup_target_name(${_targname})
  add_executable(${_targname} ${source})
  list(APPEND examples_targets ${_targname})
endforeach(source ${example_sources})

add_custom_target(examples COMMENT "Selected examples showing some of dune-gdt's features.")
add_dependencies(examples ${examples_targets})

dune_library_add_sources("boltzmann" SOURCES boltzmann-bindings.cc)
dune_library_add_sources("coordinatetransformedmn" SOURCES coordinate-transformed-mn-bindings.cc)
dune_library_add_sources("cellmodel"
                         SOURCES
                         cellmodel-bindings.cc
                         ../dune/gdt/test/cellmodel/cellmodel.cc
                         ../dune/gdt/test/cellmodel/linearsolvers.cc
                         ../dune/gdt/test/cellmodel/linear-solver-types.cc)
dune_library_add_sources("vectors" SOURCES xt-vector-bindings.cc)
