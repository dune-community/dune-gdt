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
#include "eocexpectations-fv-transport-2dyaspgrid.hh"
namespace Dune {
namespace GDT {
namespace Test {

std::vector<double> HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {7.29e-02, 7.51e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.23e-02, 1.04e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {7.34e-02, 7.64e-02, 6.26e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::godunov,
                                              TimeStepperMethods::dormand_prince>::
    results(const HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::godunov,
                                            TimeStepperMethods::dormand_prince>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {7.33e-02, 7.85e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.23e-02, 1.06e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {7.39e-02, 7.96e-02, 7.76e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                              Hyperbolic::ChooseDiscretizer::fv,
                                              2,
                                              NumericalFluxes::laxfriedrichs,
                                              TimeStepperMethods::explicit_euler>::
    results(const HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp2, double, 1>,
                                            Hyperbolic::ChooseDiscretizer::fv,
                                            2,
                                            NumericalFluxes::laxfriedrichs,
                                            TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
            const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {7.36e-02, 7.72e-02};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {1.27e-02, 1.10e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {7.41e-02, 7.85e-02, 6.67e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

} // namespace Test
} // namespace GDT
} // namespace Dune
