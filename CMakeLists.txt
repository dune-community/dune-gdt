# ~~~
# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
# License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
#      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
#          with "runtime exception" (http://www.dune-project.org/license.html)
# Authors:
#   Felix Schindler (2012 - 2018)
#   René Fritze     (2014, 2016 - 2018)
#   Sven Kaulmann   (2013)
#   Tobias Leibner  (2016, 2018)
# ~~~

# set up project
project("dune-gdt" CXX)
set(BUILD_SHARED_LIBS 1)

set(BUILD_SHARED_LIBS 1)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

cmake_minimum_required(VERSION 3.1)

if(NOT (dune-common_DIR OR dune-common_ROOT OR "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
  string(REPLACE ${CMAKE_PROJECT_NAME}
                 dune-common
                 dune-common_DIR
                 ${PROJECT_BINARY_DIR})
endif()

# dune-common
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" "${dune-common_MODULE_PATH}")
include(DuneMacros)

# dune-xt-common
find_package(dune-xt REQUIRED)
find_package(dune-pybindxi REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${dune-xt_MODULE_PATH}" "${dune-pybindxi_MODULE_PATH}")
include(DuneUtils)

# start a dune project with information from dune.module
dune_project()
if(DUNE_XT_WITH_PYTHON_BINDINGS)
  dune_enable_all_packages(INCLUDE_DIRS ${dune-gdt_SOURCE_DIR}/dune MODULE_LIBRARIES dunegdt boltzmann)
else()
  dune_enable_all_packages(INCLUDE_DIRS ${dune-gdt_SOURCE_DIR}/dune MODULE_LIBRARIES boltzmann)
endif()

add_header_listing()
make_dependent_modules_sys_included()
if(EIGEN3_FOUND)
  add_definitions("-isystem ${EIGEN3_INCLUDE_DIR}")
endif(EIGEN3_FOUND)

# search for spe10 data files
find_file(SPE10MODEL1DATA
          NAMES perm_case1.dat
          PATHS "${dune-gdt_SOURCE_DIR}/../local/src/spe10/model1"
          DOC "Location of perm_case1.dat"
          NO_DEFAULT_PATH)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)
add_subdirectory(examples EXCLUDE_FROM_ALL)

include(DunePybindxiInstallPythonPackage)
# this symlinks all files in python/ to the binary dir and install into the virtualenv from there thereby making the
# compiled extensions directly usable from the venv
dune_pybindxi_install_python_package(PATH python)
include_dependent_binary_python_dirs()
# do not change order here
add_subdirectory(python)
dxt_add_make_dependent_bindings(dune-xt)

# enable headercheck
add_definitions("-DENABLE_HEADERCHECK=1")
add_format(${CMAKE_CURRENT_SOURCE_DIR})
add_tidy(${CMAKE_CURRENT_SOURCE_DIR})
add_pylicense()

# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
