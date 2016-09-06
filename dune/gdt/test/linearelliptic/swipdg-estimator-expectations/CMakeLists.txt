# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)

# expectation targets
dune_add_library(swipdg_expectations
                NO_EXPORT
                SOURCES ao2013-2dalugrid.cxx
                        esv2007-2dalugrid.cxx
                        spe10-2dalugrid.cxx
            )

