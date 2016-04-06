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


template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {7.29e-02, 7.51e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {1.23e-02, 1.04e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {7.34e-02, 7.64e-02, 6.26e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK, anything>
    : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {7.33e-02, 7.85e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {1.23e-02, 1.04e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {7.39e-02, 7.96e-02, 7.76e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                Hyperbolic::ChooseDiscretizer::fv, 2,
                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {7.36e-02, 7.72e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0 / 5.0))
          return {1.27e-02, 1.10e-02};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {7.41e-02, 7.85e-02, 6.67e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 2,
                                         Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 2,
                                         Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

template class HyperbolicEocExpectations<Hyperbolic::TransportTestCase<Dune::YaspGrid<2>, double, 1>,
                                         Hyperbolic::ChooseDiscretizer::fv, 2,
                                         Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
