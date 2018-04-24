# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2015 - 2017)
#   Rene Milk       (2016, 2018)
#
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# ~~~

# expectation targets
dune_add_library(block_swipdg_expectations
                 SOURCES
                 block-swipdg-esv2007-2dyaspgrid.cxx)
dune_add_library(cg_expectations
                 SOURCES
                 cg-ao2013-2dyaspgrid.cxx
                 cg-er2007-2dyaspgrid.cxx
                 cg-esv2007-2dyaspgrid.cxx
                 cg-mixedboundary-2dyaspgrid.cxx
                 cg-spe10-2dyaspgrid.cxx
                 cg-ao2013-2dalugrid.cxx
                 cg-er2007-2dalugrid.cxx
                 cg-esv2007-2dalugrid.cxx
                 cg-mixedboundary-2dalugrid.cxx
                 cg-spe10-2dalugrid.cxx)
dune_add_library(swipdg_expectations
                 SOURCES
                 swipdg-ao2013-2dyaspgrid.cxx
                 swipdg-er2007-2dyaspgrid.cxx
                 swipdg-esv2007-2dyaspgrid.cxx
                 swipdg-mixedboundary-2dyaspgrid.cxx
                 swipdg-spe10-2dyaspgrid.cxx
                 swipdg-ao2013-2dalugrid.cxx
                 swipdg-er2007-2dalugrid.cxx
                 swipdg-esv2007-2dalugrid.cxx
                 swipdg-mixedboundary-2dalugrid.cxx
                 swipdg-spe10-2dalugrid.cxx)

