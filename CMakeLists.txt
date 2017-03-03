# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2012 - 2016)
#   Rene Milk       (2014)
#   Sven Kaulmann   (2013)

# set up project
project("dune-gdt" CXX)

cmake_minimum_required(VERSION 3.1)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# dune-common
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${PROJECT_SOURCE_DIR}/cmake/modules"
     "${dune-common_MODULE_PATH}")
include(DuneMacros)

find_package(dune-xt-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH
     "${dune-xt-common_MODULE_PATH}")
include(DuneUtils)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages(INCLUDE_DIRS ${dune-gdt_SOURCE_DIR}/dune
                         VERBOSE)


add_header_listing()
make_dependent_modules_sys_included()
if (EIGEN3_FOUND)
  add_definitions("-isystem ${EIGEN3_INCLUDE_DIR}")
endif (EIGEN3_FOUND)

# search for spe10 data files
find_file(SPE10MODEL1DATA
          NAMES perm_case1.dat
          PATHS "${dune-gdt_SOURCE_DIR}/../local/src/spe10/model1"
          DOC "Location of perm_case1.dat"
          NO_DEFAULT_PATH)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)

if(dune-pybindxi_FOUND)
  foreach(_file
          dune/__init__.py
          dune/gdt/__init__.py
         )
  execute_process(COMMAND ln -s ${CMAKE_CURRENT_SOURCE_DIR}/python/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})

# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
