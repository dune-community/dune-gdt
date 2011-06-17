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
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  dune/Makefile
  dune/fem/Makefile
  dune/fem/functional/Makefile
  dune/fem/functional/common/Makefile
  dune/functionals/Makefile
  dune/functionals/assembler/Makefile
  dune/functionals/assembler/local/Makefile
  dune/functionals/assembler/local/codim0/Makefile
  dune/functionals/assembler/system/Makefile
  dune/functionals/common/Makefile
  dune/functionals/constraints/Makefile
  dune/functionals/container/Makefile
  dune/functionals/discretefunctionspace/Makefile
  dune/functionals/discretefunctionspace/continuous/Makefile
  dune/functionals/discretefunctionspace/subspace/Makefile
  dune/functionals/discretefunctional/Makefile
  dune/functionals/discretefunctional/local/Makefile
  dune/functionals/discretefunctional/local/codim0/Makefile
  dune/functionals/discreteoperator/Makefile
  dune/functionals/discreteoperator/local/Makefile
  dune/functionals/discreteoperator/local/codim0/Makefile
  examples/Makefile
  examples/poisson_new/Makefile
  examples/elliptic_finite_element/Makefile
  m4/Makefile
  dune-fem-functionals.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
