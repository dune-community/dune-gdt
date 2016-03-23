# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# expectation targets
add_library(lib_test_linearelliptic_cg_discretizations_sgrid
            STATIC
            EXCLUDE_FROM_ALL
            eocexpectations-cg-ao2013-2dsgrid.cxx
            eocexpectations-cg-er2007-2dsgrid.cxx
            eocexpectations-cg-esv2007-2dsgrid.cxx
            eocexpectations-cg-mixedboundary-2dsgrid.cxx
            eocexpectations-cg-spe10-2dsgrid.cxx)
add_library(lib_test_linearelliptic_cg_discretizations_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            eocexpectations-cg-ao2013-2dalugrid.cxx
            eocexpectations-cg-er2007-2dalugrid.cxx
            eocexpectations-cg-esv2007-2dalugrid.cxx
            eocexpectations-cg-mixedboundary-2dalugrid.cxx
            eocexpectations-cg-spe10-2dalugrid.cxx)

# add some compile flags
foreach (_target lib_test_linearelliptic_cg_discretizations_alugrid)
  add_dune_alugrid_flags(${_target})
endforeach (_target)
foreach (_target lib_test_linearelliptic_cg_discretizations_alugrid
                 lib_test_linearelliptic_cg_discretizations_sgrid)
  add_dune_mpi_flags(${_target})
endforeach (_target)
