// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH

#include "discretizers/base.hh"
#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/gdt/test/grids.hh>

namespace Dune {
namespace GDT {
namespace Test {

static const auto CG = LinearElliptic::ChooseDiscretizer::cg;

namespace internal {


template <int polOrder>
class LinearEllipticEocExpectationsBase
{
public:
  static size_t rate(const std::string type)
  {
    if (type == "L2")
      return polOrder + 1;
    else if (type == "H1_semi")
      return polOrder;
    else if (type == "energy")
      return polOrder;
    else
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
    return 0;
  } // ... rate(...)
}; // class LinearEllipticEocExpectationsBase


} // namespace internal

template <class TestCaseType, LinearElliptic::ChooseDiscretizer disc, int polOrder>
class LinearEllipticEocExpectations;

} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_BASE_HH
