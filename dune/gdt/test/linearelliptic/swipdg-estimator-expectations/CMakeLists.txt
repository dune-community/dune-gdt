# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016)
# ~~~

# expectation targets
dune_add_library(swipdg_estimator_expectations
                NO_EXPORT
                SOURCES
                ao2013-2dalugrid.cxx
                esv2007-2dalugrid.cxx
                esv2007-2dyaspgrid.cxx
                spe10-2dalugrid.cxx)

