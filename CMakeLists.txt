# This file is part of the dune-gdt project:
#   https://github.com/dune-community/dune-gdt
# Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
# Authors:
#   Felix Schindler (2012 - 2016)
#   Rene Milk       (2014)
#   Sven Kaulmann   (2013)

# set up project
project("dune-gdt" C CXX)

cmake_minimum_required(VERSION 2.8.6)

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# find dune-common and set the module path
find_package(dune-common)
list(APPEND CMAKE_MODULE_PATH
     ${dune-common_MODULE_PATH}
     "${CMAKE_SOURCE_DIR}/../dune-stuff/cmake/modules"
     "${CMAKE_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)
include(AddALUGridCompileFlags)

# start a dune project with information from dune.module
dune_project()
include(DuneUtils)

file( GLOB_RECURSE gdt "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
set( COMMON_HEADER ${gdt} ${DUNE_HEADERS} )
# add header of dependent modules for header listing
foreach(_mod ${ALL_DEPENDENCIES})
file(GLOB_RECURSE HEADER_LIST "${CMAKE_CURRENT_SOURCE_DIR}/../${_mod}/*.hh")
list(APPEND COMMON_HEADER ${HEADER_LIST})
endforeach(_mod DEPENDENCIES)
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )
file( GLOB_RECURSE gdt_cxx "${CMAKE_CURRENT_SOURCE_DIR}/*.cxx" )
file( GLOB_RECURSE gdt_cc "${CMAKE_CURRENT_SOURCE_DIR}/*.cc" )
add_format(${gdt} ${gdt_cc} ${gdt_cxx})

#disable most warnings from dependent modules
foreach(_mod ${ALL_DEPENDENCIES})
dune_module_to_uppercase(_upper_case "${_mod}")
if(${_mod}_INCLUDE_DIRS)
    foreach( _idir ${${_mod}_INCLUDE_DIRS} )
    add_definitions("-isystem ${_idir}")
    endforeach( _idir )
endif(${_mod}_INCLUDE_DIRS)
endforeach(_mod DEPENDENCIES)

set( COMMON_LIBS
    ${DUNE_DEFAULT_LIBS}
    ${GRIDLIBS}
    ${FASPLIB}
)

# search for spe10 data files
find_file(SPE10MODEL1DATA
          NAMES perm_case1.dat
          PATHS "${dune-gdt_SOURCE_DIR}/../local/src/spe10/model1"
          DOC "Location of perm_case1.dat"
          NO_DEFAULT_PATH)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(cmake/modules)

# finalize the dune project, e.g., generate config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
