# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2016)
#   Tobias Leibner  (2016)

# expectation targets
dune_add_library(expectations_hyperbolic_fv_discretizations_yaspgrid
	    STATIC
            EXCLUDE_FROM_ALL
            eocexpectations-fv-sourcebeam-1dyaspgrid.cxx
            eocexpectations-fv-boltzmanncheckerboard-2dyaspgrid.cxx
            eocexpectations-fv-burgers-1dyaspgrid.cxx
            eocexpectations-fv-burgers-2dyaspgrid.cxx
            eocexpectations-fv-transport-1dyaspgrid.cxx
            eocexpectations-fv-transport-2dyaspgrid.cxx
            eocexpectations-fv-shallowwater-1dyaspgrid.cxx
            eocexpectations-fv-shocktube-1dyaspgrid.cxx)

