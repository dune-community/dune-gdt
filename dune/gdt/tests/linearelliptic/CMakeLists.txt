# This file is part of the dune-gdt project:
#   /home/falbr_01/Projects/dune/dune-hdd-demos/dune-hdd/test
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# expectation targets
add_library(obj_tests_linearelliptic_eocexpectations_cg_ao2013_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-ao2013-2dsgrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_ao2013_2dalugrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-ao2013-2dalugrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_er2007_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-er2007-2dsgrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_er2007_2dalugrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-er2007-2dalugrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_esv2007_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-esv2007-2dsgrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_esv2007_2dalugrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-esv2007-2dalugrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_mixedboundary_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-mixedboundary-2dsgrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_mixedboundary_2dalugrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-mixedboundary-2dalugrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_spe10_2dsgrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-spe10-2dsgrid.cxx)
add_library(obj_tests_linearelliptic_eocexpectations_cg_spe10_2dalugrid
            OBJECT
            EXCLUDE_FROM_ALL
            eocexpectations-cg-spe10-2dalugrid.cxx)

# test target libs
add_library(lib_test_linearelliptic_cg_discretizations_sgrid
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_ao2013_2dsgrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_er2007_2dsgrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_esv2007_2dsgrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_mixedboundary_2dsgrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_spe10_2dsgrid>)
add_library(lib_test_linearelliptic_cg_discretizations_alugrid
            STATIC
            EXCLUDE_FROM_ALL
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_ao2013_2dalugrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_er2007_2dalugrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_esv2007_2dalugrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_mixedboundary_2dalugrid>
            $<TARGET_OBJECTS:obj_tests_linearelliptic_eocexpectations_cg_spe10_2dalugrid>)
add_dune_alugrid_flags(lib_test_linearelliptic_cg_discretizations_alugrid)

