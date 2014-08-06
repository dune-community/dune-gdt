# This file is part of the dune-gdt project:
#   http://users.dune-project.org/projects/dune-gdt
# Copyright holders: Felix Schindler
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

# set up project
project("dune-gdt" C CXX)

cmake_minimum_required(VERSION 2.8.6)

# local environment
set( ENV{PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/bin:$ENV{PATH}" )
set( ENV{LD_LIBRARY_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib:$ENV{LD_LIBRARY_PATH}" )
set( ENV{PKG_CONFIG_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../local/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )


if(CMAKE_WITH_AUTOTOOLS)
  # this is how it used to be in the good old days :)
  # this assumes that we use autotools to configure everything and just use this CMakeLists.txt for the tests and
  # examples in dune-gdt
  if( NOT CMAKE_BUILD_TYPE )
    message( STATUS "No build type selected, setting default build type 'Debug'" )
    set( CMAKE_BUILD_TYPE "DEBUG" CACHE STRING "DEBUG/RELEASE")
  endif( NOT CMAKE_BUILD_TYPE )

  # cmake specific
  set( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
  set( CMAKE_COLOR_MAKEFILE ON)
  set( CMAKE_VERBOSE_MAKEFILE ON )

  # needs to come first to include correct config.h
  include_directories( ${dune-gdt_SOURCE_DIR} )
  include_directories( ${dune-gdt_BINARY_DIR} )

  # dune-stuff cmake includes
  set( CMAKE_MODULE_PATH
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake
      ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake/modules
  )
  include(DuneUtils)
  include(GridUtils)

  # we need eigen
  if( NOT EIGEN_FOUND )
    message( FATAL_ERROR "eigen3 not found, set PKG_CONFIGK_PATH to include eigen3.pc!" )
  endif( NOT EIGEN_FOUND )

  # we would like to have fasp
  #include(FindFASP)

  add_dune_modules(
    stuff
    common
    pdelab
    typetree
    fem
    fem-localfunctions
    istl
    geometry
    grid
    grid-multiscale
    localfunctions
  )

  # header
  file( GLOB_RECURSE gdt "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
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
  add_subdirectory(test)

  ADD_EXECUTABLE( this_needs_to_be_here_to_have_header_listing_until_we_have_a_real_executable "say_hello.cc" ${COMMON_HEADER} )
  TARGET_LINK_LIBRARIES( this_needs_to_be_here_to_have_header_listing_until_we_have_a_real_executable ${COMMON_LIBS} )

else(CMAKE_WITH_AUTOTOOLS)
  # this assumes that we run dunecontrol with the --use-cmake option

  # circumvent not building docs
  set(BUILD_DOCS 1)

  # find dune-common and set the module path
  find_package(dune-common)
  list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
        "${CMAKE_SOURCE_DIR}/cmake/modules")

  #include the dune macros
  include(DuneMacros)

  # start a dune project with information from dune.module
  dune_project()
  include(DuneUtils)
  set(EXTRA_DIST dune.module)
  set(SUBDIRS m4)
  if( BUILD_DOCS)
  list(APPEND SUBDIRS doc)
  endif( BUILD_DOCS )
  set(AUTOMAKE_OPTIONS  foreign 1.5)
  set(DISTCHECK_CONFIGURE_FLAGS
      --with-dune-common=${DUNE_COMMON_ROOT}
      --with-dune-grid=${DUNE_GRID_ROOT}
      --with-dune-fem=${DUNE_FEM_ROOT}
      --with-dune-stuff=${DUNE_STUFF_ROOT}
      --with-dune-geometry=${DUNE_GEOMETRY_ROOT} CXX="${CXX}" CC="${CC}")

  file( GLOB_RECURSE gdt "${CMAKE_CURRENT_SOURCE_DIR}/*.hh" )
  set( COMMON_HEADER ${gdt} ${DUNE_HEADERS} )
  set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

  set( COMMON_LIBS
      ${DUNE_DEFAULT_LIBS}
      ${GRIDLIBS}
      ${FASPLIB}
  )

  foreach(i ${SUBDIRS})
    if(${i} STREQUAL "test")
      set(opt EXCLUDE_FROM_ALL)
    endif(${i} STREQUAL "test")
    add_subdirectory(${i} ${opt})
    unset(opt)
  endforeach(i ${SUBDIRS})

  add_subdirectory(cmake/modules)
  add_subdirectory(test EXCLUDE_FROM_ALL)

  # finalize the dune project, e.g., generate config.h etc.
  finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
endif(CMAKE_WITH_AUTOTOOLS)
