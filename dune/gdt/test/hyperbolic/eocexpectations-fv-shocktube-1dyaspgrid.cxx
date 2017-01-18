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
#include "eocexpectations-fv-shocktube-1dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double> HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25))
        return {6.95e-02, 4.88e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
        return {6.68e-03, 4.69e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {6.99e-02, 4.93e-02, 3.30e-02, 2.05e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::dormand_prince>::
    results(const HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::dormand_prince>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25))
        return {7.01e-02, 4.98e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
        return {6.68e-03, 4.71e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {7.05e-02, 5.04e-02, 3.47e-02, 2.32e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              1,
                                              NumericalFluxes::laxfriedrichs,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Yasp1, double>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            1,
                                            NumericalFluxes::laxfriedrichs,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1)
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25))
        return {8.82e-02, 6.65e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
        return {8.66e-03, 6.02e-03};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    else
      return {8.86e-02, 6.70e-02, 4.91e-02, 3.42e-02};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
