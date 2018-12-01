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

#include <dune/xt/common/float_cmp.hh>

#include "eocexpectations-fv-transport-1dyaspgrid.hh"

namespace Dune {
namespace GDT {
namespace Test {


std::vector<double> Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                               Hyperbolic::ChooseDiscretizer::fv,
                                                               1,
                                                               NumericalFluxes::godunov,
                                                               TimeStepperMethods::explicit_euler,
                                                               TimeStepperMethods::explicit_euler>::
    results(
        const Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                         Hyperbolic::ChooseDiscretizer::fv,
                                                         1,
                                                         NumericalFluxes::godunov,
                                                         TimeStepperMethods::explicit_euler,
                                                         TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
        const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1)
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {3.33e-01, 2.95e-01};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {5.21e-02, 3.70e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    else
      return {3.58e-01, 3.14e-01, 2.40e-01, 1.64e-01, 1.10e-01};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                               Hyperbolic::ChooseDiscretizer::fv,
                                                               1,
                                                               NumericalFluxes::godunov,
                                                               TimeStepperMethods::dormand_prince,
                                                               TimeStepperMethods::dormand_prince>::
    results(
        const Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                         Hyperbolic::ChooseDiscretizer::fv,
                                                         1,
                                                         NumericalFluxes::godunov,
                                                         TimeStepperMethods::dormand_prince,
                                                         TimeStepperMethods::dormand_prince>::TestCaseType& test_case,
        const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1)
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {3.36e-01, 3.05e-01};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {5.23e-02, 3.74e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    else
      return {3.60e-01, 3.23e-01, 2.72e-01, 2.10e-01, 1.51e-01};
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                               Hyperbolic::ChooseDiscretizer::fv,
                                                               1,
                                                               NumericalFluxes::laxfriedrichs,
                                                               TimeStepperMethods::explicit_euler,
                                                               TimeStepperMethods::explicit_euler>::
    results(
        const Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                         Hyperbolic::ChooseDiscretizer::fv,
                                                         1,
                                                         NumericalFluxes::laxfriedrichs,
                                                         TimeStepperMethods::explicit_euler,
                                                         TimeStepperMethods::explicit_euler>::TestCaseType& test_case,
        const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {3.46e-01, 3.35e-01};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {5.91e-02, 4.77e-02};
      else
        EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
    } else {
      return {3.74e-01, 3.50e-01, 3.10e-01, 2.50e-01, 1.87e-01};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}

std::vector<double> Dune::GDT::Test::HyperbolicEocExpectations<
    Hyperbolic::TransportTestCase<Yasp1, double, 1>,
    Hyperbolic::ChooseDiscretizer::fv,
    1,
    NumericalFluxes::godunov,
    TimeStepperMethods::explicit_euler,
    TimeStepperMethods::explicit_euler,
    1>::results(const Dune::GDT::Test::HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                                                 Hyperbolic::ChooseDiscretizer::fv,
                                                                 1,
                                                                 NumericalFluxes::godunov,
                                                                 TimeStepperMethods::explicit_euler,
                                                                 TimeStepperMethods::explicit_euler,
                                                                 1>::TestCaseType& test_case,
                const std::string type)
{
  if (type == "L1") {
    if (test_case.num_refinements() == 1) {
      if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
        return {4.75e-01, 2.81e-01};
      else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
        return {4.86e-02, 2.86e-02};
      else
        EXPECT_TRUE(false) << "test results missing for type: " << type;
    } else {
      return {3.40e-01, 2.49e-01, 1.05e-01, 3.65e-02, 3.20e-02};
    }
  } else
    EXPECT_TRUE(false) << "test results missing for type: " << type;
  return {};
}


} // namespace Test
} // namespace GDT
} // namespace Dune
