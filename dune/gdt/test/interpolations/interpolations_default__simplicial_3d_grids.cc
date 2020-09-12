// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>

#include "default.hh"


using Simplicial3dGrids = ::testing::Types<
#if HAVE_DUNE_ALUGRID
    ALU_3D_SIMPLEX_CONFORMING,
    ALU_3D_SIMPLEX_NONCONFORMING
#endif
#if HAVE_DUNE_ALUGRID && HAVE_DUNE_UGGRID
    ,
#endif
#if HAVE_DUNE_UGGRID
    UG_3D
#endif
    >;


template <class G>
using InterpolationTest = Dune::GDT::Test::DefaultInterpolationOnLeafViewTest<G>;
TYPED_TEST_SUITE(InterpolationTest, Simplicial3dGrids);
TYPED_TEST(InterpolationTest, interpolates_correctly)
{
  this->interpolates_correctly(3e-14);
}
