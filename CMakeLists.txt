PROJECT( dune-detailed-discretizations)

SET( CMAKE_CXX_COMPILER
  "g++-4.6" CACHE FORCE
  "" )

# cmake settings
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0009 NEW)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
include(DuneUtils)
include(UseDoxygen)
include (FindPkgConfig)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a" )
SET( CMAKE_COLOR_MAKEFILE ON )

# user settings
SET( BLAS_LIB
  "blas" CACHE STRING
  "BLAS_LIB" )

SET( GRIDTYPE
  "SGRID" CACHE STRING
  "GRIDTYPE" )

SET( GRIDDIM
  "2" CACHE STRING
  "GRIDDIM" )

SET( POLORDER
  "1" CACHE STRING
  "POLORDER" )

SET( ALUGRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "ALUGRID_BASE_PATH" )

SET( CUSTOM_FLAGS
  "-Wall -O0 -DDEBUG -g -fdiagnostics-show-option" CACHE LIST
  "CUSTOM_FLAGS" )

# boost
set( BOOST_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../local" )
find_package( Boost 1.48.0 REQUIRED )
include_directories(${Boost_INCLUDE_DIR})
link_directories( ${Boost_LIBRARY_DIRS} )

# eigen3
set( EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../local/include/eigen3")
find_package( Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# superlu
find_package( SuperLU REQUIRED)
include_directories( ${SUPERLU_INCLUDES} )
link_directories( ${SUPERLU_LIBRARY_DIRS} )

include_dir( ${dune-detailed-discretizations_SOURCE_DIR} )
add_cxx_flags( -DHAVE_CONFIG_H )
add_cxx_flags( -D${GRIDTYPE} )
add_cxx_flags( -DGRIDDIM=${GRIDDIM} )
add_cxx_flags( -DPOLORDER=${POLORDER} )
add_cxx_flags( -I${ALUGRID_BASE_PATH}/include
  -I${ALUGRID_BASE_PATH}/include/serial
  -I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )
add_cxx_flags(-DBOOST_DATE_TIME_NO_LOCALE)
link_directories( "${ALUGRID_BASE_PATH}/lib" )
add_cxx_flags(-std=c++0x)
add_cxx_flags(${CUSTOM_FLAGS})

# other dune modules
ADD_DUNE_MODULES(common fem geometry grid grid-multiscale localfunctions stuff)

# source listing, relative to the position of this file
FILE( GLOB_RECURSE common "../dune-common/*.hh" )
FILE( GLOB_RECURSE common_cc "../dune-common/*.cc" )
FILE( GLOB_RECURSE detailed_discretizations "../dune-detailed-discretizations/*.hh" )
FILE( GLOB_RECURSE detailed_discretizations_cc "../dune-detailed-discretizations/*.cc" )
FILE( GLOB_RECURSE fem "../dune-fem/*.hh" )
FILE( GLOB_RECURSE fem_cc "../dune-fem/*.cc" )
FILE( GLOB_RECURSE geometry "../dune-geometry/*.hh" )
FILE( GLOB_RECURSE geometry_cc "../dune-geometry/*.cc" )
FILE( GLOB_RECURSE grid "../dune-grid/*.hh" )
FILE( GLOB_RECURSE grid_cc "../dune-grid/*.cc" )
FILE( GLOB_RECURSE grid_multiscale "../dune-grid-multiscale/*.hh" )
FILE( GLOB_RECURSE grid_multiscale_cc "../dune-grid-multiscale/*.cc" )
FILE( GLOB_RECURSE localfunctions "../dune-localfunctions/*.hh" )
FILE( GLOB_RECURSE localfunctions_cc "../dune-localfunctions/*.cc" )
FILE( GLOB_RECURSE stuff "../dune-stuff/*.hh" )
FILE( GLOB_RECURSE stuff_cc "../dune-stuff/*.cc" )

set_source_files_properties(
${common} ${common_cc} ${fem} ${geometry} ${grid} ${grid_multiscale} ${localfunctions} ${stuff}
${detailed_discretizations} ${detailed_discretizations_cc} ${fem_cc} ${grid_cc} ${localfunctions_cc} ${geometry_cc} ${stuff_cc}  PROPERTIES HEADER_FILE_ONLY 1 )

set(LIST_LIBS ${common} ${grid} ${grid_multiscale} ${fem} ${localfunctions} ${geometry} ${stuff} ${detailed_discretizations} )

set(COMMON_LIBS "fem" "grid" "dunegeometry" "common" ${BLAS_LIB} "alugrid" "lapack" "boost_filesystem" "boost_system" "superlu" )

ADD_EXECUTABLE( examples_elliptic_continuousgalerkin "examples/elliptic/continuousgalerkin.cc" ${LIST_LIBS} )
TARGET_LINK_LIBRARIES( examples_elliptic_continuousgalerkin ${COMMON_LIBS} )

ADD_EXECUTABLE( examples_basefunctionset_continuous_lagrange "examples/basefunctionset/continuous/lagrange.cc" ${LIST_LIBS} )
TARGET_LINK_LIBRARIES( examples_basefunctionset_continuous_lagrange ${COMMON_LIBS} )
