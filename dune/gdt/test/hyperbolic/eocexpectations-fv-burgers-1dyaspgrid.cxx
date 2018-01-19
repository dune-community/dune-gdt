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
#include "eocexpectations-fv-burgers-1dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {8.96e-02, 3.87e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.33e-02, 6.34e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {1.14e-01, 6.19e-02, 2.89e-02, 1.25e-02, 4.82e-03};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::dormand_prince,
                                              TimeStepperMethods::dormand_prince>::
    results(const HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::dormand_prince,
                                            TimeStepperMethods::dormand_prince>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {8.64e-02, 3.72e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.33e-02, 6.33e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {1.15e-01, 6.50e-02, 3.31e-02, 1.52e-02, 5.68e-03};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::laxfriedrichs,
                                              TimeStepperMethods::explicit_euler,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::BurgersTestCase<Yasp1, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::laxfriedrichs,
                                            TimeStepperMethods::explicit_euler,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1)
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {1.03e-01, 5.58e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.71e-02, 7.63e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    else
      return {1.77e-01, 1.30e-01, 7.80e-02, 3.85e-02, 1.41e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
