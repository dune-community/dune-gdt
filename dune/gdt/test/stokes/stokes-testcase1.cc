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

template <class G>
using StokesTest = StokesTestcase1<G>;

using Grids2D = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#  if HAVE_DUNE_ALUGRID
//                                 ,
//                                 ALU_2D_SIMPLEX_CONFORMING
//                                 ALU_2D_SIMPLEX_NONCONFORMING,
//                                 ALU_2D_CUBE
#  endif
#  if HAVE_DUNE_UGGRID || HAVE_UG
//                                 ,
//                                 UG_2D
#  endif
#  if HAVE_ALBERTA
//                                 ,
//                                 ALBERTA_2D
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


TYPED_TEST_CASE(StokesTest, Grids2D);

TYPED_TEST(StokesTest, order2)
{
  this->run(2);
}

#  if 0
TYPED_TEST(StokesTest, order3)
{
  this->run(3);
}
#  endif


#endif // HAVE_DUNE_ISTL
