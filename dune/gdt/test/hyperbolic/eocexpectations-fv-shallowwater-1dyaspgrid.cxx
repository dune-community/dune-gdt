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

#include "config.h"
#include "eocexpectations-fv-shallowwater-1dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::ShallowWaterTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::ShallowWaterTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 3.0))
        return {3.02e+00, 1.59e+00};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 3.0 / 5.0))
        return {5.04e-01, 2.31e-01};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {4.10e+00, 2.82e+00, 1.41e+00};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
