// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/transport.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>, Hyperbolic::ChooseDiscretizer::fv, 1,
                                NumericalFluxes::godunov, TimeStepperMethods::explicit_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::TransportTestCase<Yasp1, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
          return {3.33e-01, 2.95e-01};
        else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {5.02e-02, 3.51e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
      else
        return {3.44e-01, 3.10e-01, 2.40e-01, 1.64e-01, 1.09e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>, Hyperbolic::ChooseDiscretizer::fv, 1,
                                NumericalFluxes::godunov, TimeStepperMethods::dormand_prince, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::TransportTestCase<Yasp1, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
          return {3.36e-01, 3.05e-01};
        else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {5.03e-02, 3.55e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
      else
        return {3.46e-01, 3.19e-01, 2.72e-01, 2.10e-01, 1.51e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>, Hyperbolic::ChooseDiscretizer::fv, 1,
                                NumericalFluxes::laxfriedrichs, TimeStepperMethods::explicit_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::TransportTestCase<Yasp1, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1) {
        if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
          return {3.46e-01, 3.35e-01};
        else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {5.69e-02, 4.67e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << Dune::XT::Common::to_string(test_case.t_end());
      } else {
        return {3.57e-01, 3.48e-01, 3.12e-01, 2.50e-01, 1.87e-01};
      }
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>, Hyperbolic::ChooseDiscretizer::fv, 1,
                                NumericalFluxes::godunov_with_reconstruction, TimeStepperMethods::explicit_euler,
                                anything> : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::TransportTestCase<Yasp1, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1) {
        if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0))
          return {4.75e-01, 2.81e-01};
        else if (Dune::XT::Common::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {4.75e-02, 2.81e-02};
        else
          EXPECT_TRUE(false) << "test results missing for type: " << type;
      } else {
        return {3.29e-01, 2.47e-01, 1.06e-01, 3.83e-02, 3.33e-02};
      }
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                         TimeStepperMethods::explicit_euler>;

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::godunov,
                                         TimeStepperMethods::dormand_prince>;

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1, NumericalFluxes::laxfriedrichs,
                                         TimeStepperMethods::explicit_euler>;

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Yasp1, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1,
                                         NumericalFluxes::godunov_with_reconstruction,
                                         TimeStepperMethods::explicit_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
