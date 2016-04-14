# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

# expectation targets
add_library(expectations_linearelliptic_swipdg_estimators_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            ao2013-2dalugrid.cxx
            esv2007-2dalugrid.cxx
            spe10-2dalugrid.cxx)

# add some compile flags
foreach (_target expectations_linearelliptic_swipdg_estimators_alugrid)
  add_dune_alugrid_flags(${_target})
  add_dune_mpi_flags(${_target})
endforeach (_target)

