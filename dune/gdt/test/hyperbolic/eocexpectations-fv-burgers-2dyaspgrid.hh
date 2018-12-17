// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_BURGERS_2DYASPGRID_HH
#define DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_BURGERS_2DYASPGRID_HH

#include "problems/burgers.hh"
#include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <>
class HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp2, double, 1>,
                                Hyperbolic::ChooseDiscretizer::fv,
                                2,
                                NumericalFluxes::godunov,
                                TimeStepperMethods::explicit_euler,
                                TimeStepperMethods::explicit_euler> : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::BurgersTestCase<Yasp2, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type); // ... results(...)
}; // HyperbolicEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_BURGERS_2DYASPGRID_HH
