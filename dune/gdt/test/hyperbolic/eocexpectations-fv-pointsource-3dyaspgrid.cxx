// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#include <config.h>
#include "eocexpectations-fv-pointsource-3dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::PointSourceTestCase<Yasp3, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              3,
                                              NumericalFluxes::laxfriedrichs,
                                              TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                              TimeStepperMethods::matrix_exponential>::
    results(const HyperbolicEocExpectations<Hyperbolic::PointSourceTestCase<Yasp3, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            3,
                                            NumericalFluxes::laxfriedrichs,
                                            TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                            TimeStepperMethods::matrix_exponential>::TestCaseType& /*test_case*/,
            const std::string type)
{
  if (type == "L1") {
    return {3.25e-01, 1.63e-01};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
