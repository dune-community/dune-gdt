// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#define DUNE_XT_COMMON_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/stokes/stokes-taylorhood.hh>

#if HAVE_DUNE_ISTL

using namespace Dune;
using namespace Dune::GDT::Test;

using SimplexGrids2D = ::testing::Types<
#  if HAVE_DUNE_ALUGRID
    ALU_2D_SIMPLEX_CONFORMING,
    ALU_2D_SIMPLEX_NONCONFORMING
#  endif
#  if HAVE_DUNE_UGGRID || HAVE_UG
#    if HAVE_DUNE_ALUGRID
    ,
#    endif
    UG_2D
#  endif
    >;

using CubeGrids2D = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#  if HAVE_DUNE_ALUGRID
                                     ,
                                     ALU_2D_CUBE
#  endif
                                     >;

DUNE_XT_COMMON_TYPENAME(YASP_2D_EQUIDISTANT_OFFSET)
#  if HAVE_DUNE_ALUGRID
DUNE_XT_COMMON_TYPENAME(ALU_2D_SIMPLEX_CONFORMING)
DUNE_XT_COMMON_TYPENAME(ALU_2D_SIMPLEX_NONCONFORMING)
DUNE_XT_COMMON_TYPENAME(ALU_2D_CUBE)
#  endif
#  if HAVE_DUNE_UGGRID || HAVE_UG
DUNE_XT_COMMON_TYPENAME(UG_2D)
#  endif
#  if HAVE_ALBERTA
DUNE_XT_COMMON_TYPENAME(ALBERTA_2D)
#  endif


template <class G>
using StokesTestSimplex = StokesTestcase1<G>;
TYPED_TEST_CASE(StokesTestSimplex, SimplexGrids2D);

TYPED_TEST(StokesTestSimplex, order2)
{
  this->run(2, 2e-5, 3e-3);
}

template <class G>
using StokesTestCube = StokesTestcase1<G>;
TYPED_TEST_CASE(StokesTestCube, CubeGrids2D);

TYPED_TEST(StokesTestCube, order2)
{
  this->run(2, 3e-6, 3e-5);
}

TYPED_TEST(StokesTestCube, order3)
{
  this->run(3, 4e-7, 7e-6);
}

#endif // HAVE_DUNE_ISTL
