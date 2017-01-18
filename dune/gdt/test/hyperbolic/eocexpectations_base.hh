// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH
#define DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH

#include <dune/xt/common/test/gtest/gtest.h>

#include "discretizers/base.hh"

#include <dune/gdt/operators/fv.hh>
#include <dune/gdt/test/grids.hh>
#include <dune/gdt/timestepper/interface.hh>

namespace Dune {
namespace GDT {
namespace Test {


using Yasp1 = Yasp1Grid;
using Yasp2 = Yasp2Grid;

namespace internal {


template <int dimDomain>
class HyperbolicEocExpectationsBase
{
public:
  static double rate(const std::string type)
  {
    if (type == "L1") {
      if (dimDomain == 1)
        return 0.5;
      else
        return 0.25;
    } else {
      EXPECT_TRUE(false) << "expected rate missing for type: " << type;
      return 0;
    }
  } // ... rate(...)
}; // class HyperbolicEocExpectationsBase


} // namespace internal

template <class TestCaseType,
          Hyperbolic::ChooseDiscretizer disc,
          size_t dimDomain,
          NumericalFluxes num_flux,
          TimeStepperMethods time_stepper>
class HyperbolicEocExpectations;


} // namespace Test
} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_BASE_HH
