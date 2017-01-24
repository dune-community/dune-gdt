# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2010, 2013 - 2016)

if(dune-pybindxi_FOUND)
  dune_pybindxi_add_module(_gdt bindings.cc)
  target_link_dune_default_libraries(_gdt)
endif()

add_subdirectory(test EXCLUDE_FROM_ALL)

