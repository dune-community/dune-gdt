// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#include <config.h>
#include "eocexpectations-fv-burgers-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double> HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp2, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp2, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
      return {1.79e-03, 8.44e-04};
    else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
      return {3.60e-04, 1.70e-04};
    else
      EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
  } else {
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  }
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
