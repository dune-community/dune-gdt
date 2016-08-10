# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2015 - 2016)

# expectation targets
# CG
add_library(expectations_linearelliptic_cg_discretizations_yaspgrid
            STATIC
            EXCLUDE_FROM_ALL
            cg-ao2013-2dyaspgrid.cxx
            cg-er2007-2dyaspgrid.cxx
            cg-esv2007-2dyaspgrid.cxx
            cg-mixedboundary-2dyaspgrid.cxx
            cg-spe10-2dyaspgrid.cxx)
add_library(expectations_linearelliptic_cg_discretizations_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            cg-ao2013-2dalugrid.cxx
            cg-er2007-2dalugrid.cxx
            cg-esv2007-2dalugrid.cxx
            cg-mixedboundary-2dalugrid.cxx
            cg-spe10-2dalugrid.cxx)

# SWIPDG
add_library(expectations_linearelliptic_swipdg_discretizations_yaspgrid
            STATIC
            EXCLUDE_FROM_ALL
            swipdg-ao2013-2dyaspgrid.cxx
            swipdg-er2007-2dyaspgrid.cxx
            swipdg-esv2007-2dyaspgrid.cxx
            swipdg-mixedboundary-2dyaspgrid.cxx
            swipdg-spe10-2dyaspgrid.cxx)
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
                 expectations_linearelliptic_cg_discretizations_yaspgrid
                 expectations_linearelliptic_swipdg_discretizations_alugrid
                 expectations_linearelliptic_swipdg_discretizations_yaspgrid)
  add_dune_mpi_flags(${_target})
endforeach (_target)
