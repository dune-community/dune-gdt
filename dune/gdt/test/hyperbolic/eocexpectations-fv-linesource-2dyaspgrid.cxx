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
#include "eocexpectations-fv-linesource-2dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::LineSourceTestCase<Yasp2, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::laxfriedrichs,
                                              TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                              TimeStepperMethods::matrix_exponential>::
    results(const HyperbolicEocExpectations<Hyperbolic::LineSourceTestCase<Yasp2, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::laxfriedrichs,
                                            TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                            TimeStepperMethods::matrix_exponential>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.05))
      return {3.25e-01, 1.63e-01, 1.3e-01, 1.2e-01, 1.1e-01, 1.2e-01, 1.1e-01};
    else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 4.0 / 5.0))
      return {9.16e-02, 4.08e-02};
    else
      EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
