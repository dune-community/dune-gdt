// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/transport.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperKombinations::godunov_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        return {3.33e-01, 2.95e-01};
      else
        return {3.44e-01, 3.10e-01, 2.40e-01, 1.64e-01, 1.09e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperKombinations::godunov_adaptiveRK,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        return {3.36e-01, 3.05e-01};
      else
        return {3.46e-01, 3.19e-01, 2.72e-01, 2.10e-01, 1.51e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperKombinations::laxfriedrichs_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        return {3.46e-01, 3.35e-01};
      else
        return {3.57e-01, 3.48e-01, 3.12e-01, 2.50e-01, 1.87e-01};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperKombinations::godunovwithreconstruction_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        return {3.18e-01, 2.32e-01};
      else
        return {3.29e-01, 2.47e-01, 1.06e-01, 3.83e-02, 3.33e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperKombinations::godunov_euler >;

template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperKombinations::godunov_adaptiveRK >;

template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperKombinations::laxfriedrichs_euler >;

template class HyperbolicEocExpectations< Hyperbolic::TransportTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperKombinations::godunovwithreconstruction_euler >;


} // namespace Tests
} // namespace GDT
} // namespace Dune
