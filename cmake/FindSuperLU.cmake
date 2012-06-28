if (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)
  set(SUPERLU_FIND_QUIETLY TRUE)
endif (SUPERLU_INCLUDES AND SUPERLU_LIBRARIES)

find_path(SUPERLU_INCLUDES
  NAMES supermatrix.h
  HINTS ${SUPERLU_ROOT}
  PATH_SUFFIXES "SRC" "include"
  )

find_path(SUPERLU_LIBRARY_DIRS
  NAMES "libsuperlu.a" "libsuperlu_4.3.a"
  HINTS ${SUPERLU_ROOT}
  PATH_SUFFIXES "lib"
  )
  
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUPERLU DEFAULT_MSG
                                  SUPERLU_INCLUDES SUPERLU_LIBRARY_DIRS)

mark_as_advanced(SUPERLU_INCLUDES SUPERLU_LIBRARY_DIRS)
