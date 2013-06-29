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

# needs to come first to include correct config.h
include_dir( ${dune-gdt_SOURCE_DIR} )
include_dir( ${dune-gdt_BINARY_DIR} )

# user settings
set( ENABLE_MPI
     "0" CACHE STRING
     "set to 1 to enable mpi")

set( CUSTOM_FLAGS
     "-Wall -Wextra -O0 -DDEBUG -g -ggdb -fdiagnostics-show-option" CACHE LIST
     "CUSTOM_FLAGS" )

set( POLORDER
     "1" CACHE STRING
     "POLORDER" )

# we need eigen
if( NOT EIGEN_FOUND )
  message( FATAL_ERROR "eigen3 not found, set PKG_CONFIGK_PATH to include eigen3.pc!" )
endif( NOT EIGEN_FOUND )

add_definitions( ${CUSTOM_FLAGS} )

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

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
     ${PARALIBS}
     ${CCGNU_LIBRARIES}
     "boost_system" "boost_filesystem"
     ${GRIDLIBS}
     ${FASPLIB}
)

add_subdirectory(doc)
add_subdirectory(dune)
add_subdirectory(examples)

