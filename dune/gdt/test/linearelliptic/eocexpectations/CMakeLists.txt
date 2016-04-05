# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# expectation targets
# CG
add_library(expectations_linearelliptic_cg_discretizations_sgrid
            STATIC
            EXCLUDE_FROM_ALL
            cg-ao2013-2dsgrid.cxx
            cg-er2007-2dsgrid.cxx
            cg-esv2007-2dsgrid.cxx
            cg-mixedboundary-2dsgrid.cxx
            cg-spe10-2dsgrid.cxx)
add_library(expectations_linearelliptic_cg_discretizations_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            cg-ao2013-2dalugrid.cxx
            cg-er2007-2dalugrid.cxx
            cg-esv2007-2dalugrid.cxx
            cg-mixedboundary-2dalugrid.cxx
            cg-spe10-2dalugrid.cxx)

# SWIPDG
add_library(expectations_linearelliptic_swipdg_discretizations_sgrid
            STATIC
            EXCLUDE_FROM_ALL
            swipdg-ao2013-2dsgrid.cxx
            swipdg-er2007-2dsgrid.cxx
            swipdg-esv2007-2dsgrid.cxx
            swipdg-mixedboundary-2dsgrid.cxx
            swipdg-spe10-2dsgrid.cxx)
add_library(expectations_linearelliptic_swipdg_discretizations_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            swipdg-ao2013-2dalugrid.cxx
            swipdg-er2007-2dalugrid.cxx
            swipdg-esv2007-2dalugrid.cxx
            swipdg-mixedboundary-2dalugrid.cxx
            swipdg-spe10-2dalugrid.cxx)

# add some compile flags
foreach (_target expectations_linearelliptic_cg_discretizations_alugrid
                 expectations_linearelliptic_swipdg_discretizations_alugrid)
  add_dune_alugrid_flags(${_target})
endforeach (_target)
foreach (_target expectations_linearelliptic_cg_discretizations_alugrid
                 expectations_linearelliptic_cg_discretizations_sgrid
                 expectations_linearelliptic_swipdg_discretizations_alugrid
                 expectations_linearelliptic_swipdg_discretizations_sgrid)
  add_dune_mpi_flags(${_target})
endforeach (_target)
