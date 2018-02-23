# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2016 - 2017)
#   Rene Milk       (2016, 2018)
#   Tobias Leibner  (2016 - 2017)

# expectation targets
dune_add_library(expectations_hyperbolic_fv_discretizations_yaspgrid
                 SOURCES
                 eocexpectations-fv-sourcebeam-1dyaspgrid.cxx
                 eocexpectations-fv-boltzmanncheckerboard-2dyaspgrid.cxx
                 eocexpectations-fv-burgers-1dyaspgrid.cxx
                 eocexpectations-fv-burgers-2dyaspgrid.cxx
                 eocexpectations-fv-transport-1dyaspgrid.cxx
                 eocexpectations-fv-transport-2dyaspgrid.cxx
                 eocexpectations-fv-shallowwater-1dyaspgrid.cxx
                 eocexpectations-fv-shocktube-1dyaspgrid.cxx
)
