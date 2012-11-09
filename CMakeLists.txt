PROJECT( dune-detailed-discretizations)

SET( CMAKE_CXX_COMPILER
  "g++-4.6" CACHE FORCE
  "" )

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)
cmake_policy(SET CMP0009 NEW)
# tell cmake where to look for .cmake files
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_CURRENT_SOURCE_DIR}/../dune-stuff/cmake )

# stuff needed for dune
include(DuneUtils)

# General Settings
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

SET( ALBERTA_GRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "ALBERTA_GRID_BASE_PATH" )

SET( UGGRID_BASE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../local" CACHE STRING
  "UGGRID_BASE_PATH" )

SET( CUSTOM_FLAGS
  "-Wall -O0 -DDEBUG -g -fdiagnostics-show-option" CACHE LIST
  "CUSTOM_FLAGS" )

# superlu
find_package( SuperLU REQUIRED)
include_directories( ${SUPERLU_INCLUDES} )
link_directories( ${SUPERLU_LIBRARY_DIRS} )

# eigen3
set( EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../local/include/eigen3")
find_package( Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

SET_CONFIGHEADER_VARS()
# needs to come first to include correct config.h
include_dir( ${dune-detailed-discretizations_SOURCE_DIR} )
INCLUDE_DIR( ${dune-detailed-discretizations_BINARY_DIR} )

IF( ${GRIDTYPE} STREQUAL "UGGRID" )
  ADD_CXX_FLAGS( -DENABLE_UG )
  INCLUDE_SYS_DIR( ${UGGRID_BASE_PATH}/include )
  LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
  SET( GRIDLIBS "ugS3" "ugS2" "devS")
ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

SET( ENABLE_ALUGRID OFF )
IF( ${GRIDTYPE} MATCHES "ALUGRID" )
  INCLUDE_SYS_DIR(  ${GRAPE_PATH} 
    ${ALUGRID_BASE_PATH}/include 
    ${ALUGRID_BASE_PATH}/include/serial
    ${ALUGRID_BASE_PATH}/include/duneinterface )
  LINK_DIRECTORIES( "${ALUGRID_BASE_PATH}/lib" )
  SET( GRIDLIBS "alugrid" )
  SET( ENABLE_ALUGRID ON )
  IF ( ENABLE_PARALELL )
    INCLUDE_SYS_DIR( ${ALUGRID_BASE_PATH}/include/parallel )
  ENDIF ( ENABLE_PARALELL )
  ADD_CXX_FLAGS( "-DENABLE_ALUGRID" )
ENDIF( ${GRIDTYPE} MATCHES "ALUGRID" )

SET( ENABLE_ALBERTA OFF )
IF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )
  INCLUDE_SYS_DIR(  ${ALBERTA_GRID_BASE_PATH}/include 
    ${ALBERTA_GRID_BASE_PATH}/include/alberta )
  LINK_DIRECTORIES( "${ALBERTA_GRID_BASE_PATH}/lib" )
  SET( GRIDLIBS alberta_utilities alberta_2d dunealbertagrid_2d dunegrid)
  SET( ENABLE_ALBERTA ON )
  ADD_CXX_FLAGS( -DENABLE_ALBERTA=1 -DDEBUG )
ENDIF( ${GRIDTYPE} MATCHES "ALBERTAGRID" )

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/cmake_config.h )

# other dune modules
ADD_DUNE_MODULES(common fem geometry grid grid-multiscale localfunctions stuff )
set( COMMON_LIBS "dunefem" "dunegrid" "dunecommon" "dunegeometry"
  ${BLAS_LIB} ${GRIDLIBS} ${ParaLIBS} 
  "boost_filesystem" "boost_system" )
FILE( GLOB_RECURSE header "${CMAKE_CURRENT_SOURCE_DIR}/dune/*.hh" )
set( COMMON_HEADER ${header} ${DUNE_HEADERS} )

set_source_files_properties( ${COMMON_HEADER} PROPERTIES HEADER_FILE_ONLY 1 )

add_subdirectory(doc)
add_subdirectory(examples)

HEADERCHECK( ${header} )
DEPENDENCYCHECK( ${header} )

