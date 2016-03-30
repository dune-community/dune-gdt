# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# expectation targets
add_library(expectations_hyperbolic_fv_discretizations_yaspgrid
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

# add some compile flags
foreach (_target expectations_hyperbolic_fv_discretizations_yaspgrid)
  add_dune_mpi_flags(${_target})
endforeach (_target)

