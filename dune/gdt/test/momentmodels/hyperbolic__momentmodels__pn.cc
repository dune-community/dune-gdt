// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Tobias Leibner  (2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/test/main.hxx>

#define USE_FULL_LINEAR_RECONSTRUCTION_OPERATOR 1
#include <dune/gdt/test/momentmodels/hyperbolic_momentmodels_pn_base.hh>
#include <dune/gdt/test/momentmodels/pn-discretization.hh>

TYPED_TEST_CASE(HyperbolicPnTest, YaspGridTestCasesWithoutReconstruction);
TYPED_TEST(HyperbolicPnTest, check)
{
  this->run();
}

template <class TestCaseType>
struct HyperbolicPnTestWithReconstruction : public HyperbolicPnTest<TestCaseType>
{};

TYPED_TEST_CASE(HyperbolicPnTestWithReconstruction, YaspGridTestCasesWithReconstruction);
TYPED_TEST(HyperbolicPnTestWithReconstruction, check_with_full_linear_reconstruction)
{
  // evaluation of reconstructed function is not thread-safe at the moment
  DXTC_CONFIG["threading.max_count"] = "1";
  this->run();
}
