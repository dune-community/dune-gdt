# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([])
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
  dune/detailed-discretizations/Makefile
  dune/detailed-discretizations/assembler/Makefile
  dune/detailed-discretizations/assembler/local/Makefile
  dune/detailed-discretizations/assembler/local/codim0/Makefile
  dune/detailed-discretizations/assembler/local/codim1/Makefile
  dune/detailed-discretizations/assembler/local/combined/Makefile
  dune/detailed-discretizations/assembler/system/Makefile
  dune/detailed-discretizations/basefunctionset/Makefile
  dune/detailed-discretizations/basefunctionset/continuous/Makefile
  dune/detailed-discretizations/basefunctionset/local/Makefile
  dune/detailed-discretizations/constraints/Makefile
  dune/detailed-discretizations/container/Makefile
  dune/detailed-discretizations/container/factory/Makefile
  dune/detailed-discretizations/discretefunction/Makefile
  dune/detailed-discretizations/discretefunctional/Makefile
  dune/detailed-discretizations/discretefunctional/local/Makefile
  dune/detailed-discretizations/discretefunctional/local/codim0/Makefile
  dune/detailed-discretizations/discretefunctional/local/codim1/Makefile
  dune/detailed-discretizations/discretefunctionspace/Makefile
  dune/detailed-discretizations/discretefunctionspace/continuous/Makefile
  dune/detailed-discretizations/discretefunctionspace/discontinuous/Makefile
  dune/detailed-discretizations/discretefunctionspace/subspace/Makefile
  dune/detailed-discretizations/discreteoperator/Makefile
  dune/detailed-discretizations/discreteoperator/local/Makefile
  dune/detailed-discretizations/discreteoperator/local/codim0/Makefile
  dune/detailed-discretizations/discreteoperator/local/codim1/Makefile
  dune/detailed-discretizations/evaluation/Makefile
  dune/detailed-discretizations/evaluation/local/Makefile
  dune/detailed-discretizations/evaluation/local/binary/Makefile
  dune/detailed-discretizations/evaluation/local/quaternary/Makefile
  dune/detailed-discretizations/evaluation/local/unary/Makefile
  dune/detailed-discretizations/mapper/Makefile
  dune/detailed-discretizations/mapper/continuous/Makefile
  examples/Makefile
  examples/elliptic_finite_element/Makefile
  examples/elliptic_discontinuous_galerkin/Makefile
  m4/Makefile
  dune-detailed-discretizations.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
