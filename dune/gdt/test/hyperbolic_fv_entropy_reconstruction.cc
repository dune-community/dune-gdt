// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2019)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include <dune/gdt/test/fv-discretization.hh>
#include <dune/gdt/test/momentmodels/isentropic_euler.hh>

using Yasp1 = Dune::YaspGrid<1, Dune::EquidistantOffsetCoordinates<double, 1>>;

GTEST_TEST(FvEntropyTest, check)
{
  using Problem = Dune::GDT::IsentropicEulerEquations<typename Yasp1::template Codim<0>::Entity>;
  HyperbolicFvDiscretization<Yasp1, Problem, true> test;
  test.run();
}
