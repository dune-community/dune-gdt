PROJECT(dune-fem-functionals)

#set minimum cmake version
cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
include(DuneUtils)
include(UseDoxygen)

# FIND_PACKAGE(Boost 1.35.0 REQUIRED)
# INCLUDE_DIR(${Boost_INCLUDE_DIR})
# FIND_PACKAGE( PkgConfig )
# pkg_check_modules( CCGNU REQUIRED libccgnu2 )
# ADD_CXX_FLAGS( "${CCGNU_CFLAGS}" )

#----------------------------------------------------------------------------------------------------
# general Settings
#----------------------------------------------------------------------------------------------------
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true )
SET( CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".lib" ".la" ".a")
SET( CMAKE_COLOR_MAKEFILE ON)

# EXECUTE_PROCESS(
# 	COMMAND ${CMAKE_SOURCE_DIR}/getRevision.sh
# 	OUTPUT_VARIABLE COMMIT
# 	ERROR_VARIABLE shell_error
# 	OUTPUT_STRIP_TRAILING_WHITESPACE
# )


#----------------------------------------------------------------------------------------------------
# options, that can be changed by the user in order to customise the dune module
#----------------------------------------------------------------------------------------------------
SET( BLAS_LIB
  "blas" CACHE STRING
  "BLAS_LIB" )

SET( GRIDTYPE
  "YASPGRID" CACHE STRING
  "GRIDTYPE" )

SET( GRIDDIM
  "2" CACHE STRING
  "GRIDDIM" )

SET( POLORDER
  "1" CACHE STRING
  "POLORDER" )

SET( ALUGRID_BASE_PATH
  "/opt/dune/modules/alugrid/alugrid-1.22/" CACHE STRING
  "ALUGRID_BASE_PATH" )
#
# SET( UGGRID_BASE_PATH
# 	"/share/dune/Modules/modules_x86_64/ug" CACHE STRING
# 	"UGGRID_BASE_PATH" )
#
# SET( GRAPE_PATH
# 	"/share/dune/Modules/modules_x86_64/grape" CACHE STRING
# 	"GRAPE_PATH" )

# SET( CXX_FLAGS
#   -g -Wall -Wfatal-errors CACHE STRING
#   "CXX_FLAGS")

SET( CUSTOM_FLAGS
	"-Wall -O0 -DDEBUG -funroll-loops -g -ggdb -fno-strict-aliasing -std=c++0x -fdiagnostics-show-option " CACHE LIST
	"CUSTOM FLAGS" )

# SET( ENABLE_PARALELL
# 	0 CACHE BOOL
# 	"Enable parallel features" )

# SET ( METIS_DIR
# 	"/share/dune/Modules/modules_x86_64/ParMetis-3.1.1" CACHE STRING
# 	"metis toplevel directory" )

#needs to come first to include correct config.h
INCLUDE_DIR( ${dune-fem-functionals_SOURCE_DIR} )

# IF( ${GRIDTYPE} STREQUAL "UGGRID" )
# 	ADD_CXX_FLAGS( -I${UGGRID_BASE_PATH}/include -DENABLE_UG )
# 	LINK_DIRECTORIES( ${UGGRID_BASE_PATH}/lib/ )
# 	SET( UGLIBS "ugS3" "ugS2" "devS")
# ENDIF( ${GRIDTYPE} STREQUAL "UGGRID" )

ADD_CXX_FLAGS(-DHAVE_CONFIG_H )
ADD_CXX_FLAGS(-D${GRIDTYPE} -DGRIDDIM=${GRIDDIM} )
ADD_CXX_FLAGS(-DPOLORDER=${POLORDER} )
ADD_CXX_FLAGS( ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS} )
ADD_CXX_FLAGS( -I${ALUGRID_BASE_PATH}/include
  -I${ALUGRID_BASE_PATH}/include/serial
  -I${ALUGRID_BASE_PATH}/include/duneinterface -DENABLE_ALUGRID )

# IF ( ENABLE_PARALELL )
# 	LINK_DIRECTORIES( ${METIS_DIR} )
# 	Set ( ParaLIBS "-lparmetis -lmetis" )
# 	ADD_CXX_FLAGS( -I${ALUGRID_BASE_PATH}/include/parallel
# 					-DMPIPP_H -DENABLE_MPI=1
# 					)
#
# ENDIF ( ENABLE_PARALELL )

ADD_DUNE_MODULES( common istl grid fem fem-tools)

LINK_DIRECTORIES( "/usr/lib" "${ALUGRID_BASE_PATH}/lib" )
ADD_CXX_FLAGS( ${CUSTOM_FLAGS} )

#----------------------------------------------------------------------------------------------------
# Source listing
#----------------------------------------------------------------------------------------------------
FILE( GLOB_RECURSE common "../dune-common/dune/common/*.hh" )
FILE( GLOB_RECURSE grid "../dune-grid/dune/grid/*.hh" )
FILE( GLOB_RECURSE istl "../dune-istl/dune/istl/*.hh" )
FILE( GLOB_RECURSE femhowto "../dune-fem-howto/dune/fem-howto/*.hh" )
FILE( GLOB_RECURSE fem "../dune-fem/dune/fem/*.hh" )
FILE( GLOB_RECURSE femtools "../dune-fem-tools/dune/fem-tools/*.hh" )
FILE( GLOB_RECURSE femfunctionals "../dune-fem-functionals/dune/fem/*.hh" )

set_source_files_properties( ${common} ${grid} ${istl} ${fem} ${femhowto} ${femtools} ${femfunctionals} PROPERTIES HEADER_FILE_ONLY 1 )
set( COMMON_LIBS "fem" "grid" "common" ${BLAS_LIB} "alugrid" "GL" )

ADD_EXECUTABLE( ltwofunctional_test "tests/functional/ltwofunctional_test.cc" ${common} ${grid} ${istl} ${fem} ${femhowto} ${femtools} ${femfunctionals} )
TARGET_LINK_LIBRARIES( ltwofunctional_test ${COMMON_LIBS}  )

ADD_EXECUTABLE( ellipticfiniteelementoperator_test "tests/operator/ellipticfiniteelementoperator_test.cc" ${common} ${grid} ${istl} ${fem} ${femhowto} ${femtools} ${femfunctionals} )
TARGET_LINK_LIBRARIES( ellipticfiniteelementoperator_test ${COMMON_LIBS}  )

ADD_EXECUTABLE( localbasefunction_test "tests/common/localbasefunction_test.cc" ${common} ${grid} ${istl} ${fem} ${femhowto} ${femtools} ${femfunctionals} )
TARGET_LINK_LIBRARIES( localbasefunction_test ${COMMON_LIBS}  )

ADD_EXECUTABLE( ellipticfiniteelement_example "examples/finite_element/main.cc" ${common} ${grid} ${istl} ${fem} ${femhowto} ${femtools} ${femfunctionals} )
TARGET_LINK_LIBRARIES( ellipticfiniteelement_example ${COMMON_LIBS}  )

#HEADERCHECK( ${header} )
