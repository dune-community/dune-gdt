# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

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

