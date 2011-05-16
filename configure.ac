# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/dune_fem_functionals.cc])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-fem dune-istl
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-fem]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/fem/Makefile
  dune/fem/common/Makefile
  dune/fem/constraints/Makefile
  dune/fem/container/Makefile
  dune/fem/functional/common/Makefile
  dune/fem/functional/Makefile
  dune/fem/operator/Makefile
  dune/fem/solver/Makefile
  dune/fem/subspace/Makefile
  examples/Makefile
  examples/poisson_new/Makefile
  examples/finite_element_local_operation/Makefile
  tests/Makefile
  tests/common/Makefile
  m4/Makefile
  dune-fem-functionals.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
