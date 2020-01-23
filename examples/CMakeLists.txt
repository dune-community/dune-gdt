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

macro(cleanup_target_name arg)
  string(REPLACE ${PROJECT_SOURCE_DIR}
                 ""
                 targname
                 ${arg})
  string(REPLACE ".cc"
                 ""
                 targname
                 ${targname})
  string(REGEX
         REPLACE "/"
                 "__"
                 targname
                 ${targname})
  string(REGEX
         REPLACE "^__"
                 ""
                 targname
                 ${targname})
endmacro(cleanup_target_name)

file(GLOB_RECURSE example_sources "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
foreach(source ${example_sources}) # do some name conversion
  set(targname ${source})
  cleanup_target_name(${targname})
  add_executable(${targname} ${source})
  list(APPEND examples_targets ${targname})
endforeach(source ${example_sources})

add_custom_target(examples)
add_dependencies(examples ${examples_targets})

dune_library_add_sources("boltzmann" SOURCES boltzmann-bindings.cc)
dune_library_add_sources("cellmodel" SOURCES cellmodel-bindings.cc ../dune/gdt/test/cellmodel/cellmodel.cc)
dune_library_add_sources("vectors" SOURCES xt-vector-bindings.cc)
