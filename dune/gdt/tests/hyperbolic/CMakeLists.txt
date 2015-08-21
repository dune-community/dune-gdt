# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# expectation targets
add_library(obj_tests_hyperbolic_eocexpectations_fv_burgers_1dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-burgers-1dyaspgrid.cxx)
add_library(obj_tests_hyperbolic_eocexpectations_fv_burgers_2dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-burgers-2dyaspgrid.cxx)
add_library(obj_tests_hyperbolic_eocexpectations_fv_transport_1dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-transport-1dyaspgrid.cxx)
add_library(obj_tests_hyperbolic_eocexpectations_fv_transport_2dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-transport-2dyaspgrid.cxx)
add_library(obj_tests_hyperbolic_eocexpectations_fv_shallowwater_1dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-shallowwater-1dyaspgrid.cxx)
add_library(obj_tests_hyperbolic_eocexpectations_fv_shocktube_1dyaspgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-fv-shocktube-1dyaspgrid.cxx)
# add alugrid compile flags for the above objects (add_dune_alugrid_flags does not work here)
foreach(_target obj_tests_hyperbolic_eocexpectations_fv_burgers_1dyaspgrid
		obj_tests_hyperbolic_eocexpectations_fv_burgers_2dyaspgrid
		obj_tests_hyperbolic_eocexpectations_fv_transport_1dyaspgrid
		obj_tests_hyperbolic_eocexpectations_fv_transport_2dyaspgrid
		obj_tests_hyperbolic_eocexpectations_fv_shallowwater_1dyaspgrid
		obj_tests_hyperbolic_eocexpectations_fv_shocktube_1dyaspgrid)
  add_dune_alugrid_compile_flags(${_target})
endforeach(_target)

# test target libs
add_library(lib_test_hyperbolic_fv_discretizations_yaspgrid
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_burgers_1dyaspgrid>
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_burgers_2dyaspgrid>
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_transport_1dyaspgrid>
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_transport_2dyaspgrid>
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_shallowwater_1dyaspgrid>
            $<TARGET_OBJECTS:obj_tests_hyperbolic_eocexpectations_fv_shocktube_1dyaspgrid>)

