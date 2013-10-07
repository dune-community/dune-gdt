# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Albrecht
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

if( NOT CMAKE_BUILD_TYPE )
  message( STATUS "No build type selected, setting default build type 'Debug'" )
  set( CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "DEBUG/RELEASE")
endif( NOT CMAKE_BUILD_TYPE )

project( dune-gdt )

# cmake specific
cmake_minimum_required( VERSION 2.8 )
set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
set( CMAKE_COLOR_MAKEFILE ON)
set( CMAKE_VERBOSE_MAKEFILE ON )

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )

# dune-stuff cmake includes
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )
include(DuneUtils)
include(GridUtils)

# adjust CXX_FLAGS
set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3" )
set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -g3 -ggdb" )

# needs to come first to include correct config.h
include_dir( ${dune-gdt_SOURCE_DIR} )
include_dir( ${dune-gdt_BINARY_DIR} )

# we need eigen
if( NOT EIGEN_FOUND )
  message( FATAL_ERROR "eigen3 not found, set PKG_CONFIGK_PATH to include eigen3.pc!" )
endif( NOT EIGEN_FOUND )

# we would like to have fasp
include(FindFASP)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )
add_definitions( -DHAVE_CMAKE_CONFIG )

add_dune_modules(
  stuff
  common
  fem
  fem-localfunctions
  geometry
  grid
  grid-multiscale
  localfunctions
)

# header
file( GLOB_RECURSE gdt "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${gdt} ${DUNE_HEADERS} )
set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

# libs
set( COMMON_LIBS
     ${DUNE_LIBS}
     ${DUNE_DEFAULT_LIBS}
     ${PARALIBS}
     ${CCGNU_LIBRARIES}
     ${BOOST_SYSTEM_LIBS}
     ${GRIDLIBS}
     ${FASPLIB}
)

add_subdirectory(doc)
add_subdirectory(dune)

