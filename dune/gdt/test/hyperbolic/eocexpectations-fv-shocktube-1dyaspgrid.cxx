// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/sodshocktube.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                Hyperbolic::FluxTimeStepperCombinations::godunov_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 0.25))
          return {6.95e-02, 4.88e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
          return {6.68e-03, 4.69e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {6.99e-02, 4.93e-02, 3.30e-02, 2.05e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 0.25))
          return {7.01e-02, 4.98e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
          return {6.68e-03, 4.69e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {7.05e-02, 5.04e-02, 3.47e-02, 2.32e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template <bool anything>
class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                Hyperbolic::ChooseDiscretizer::fv, 1,
                                Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler, anything>
    : public internal::HyperbolicEocExpectationsBase<1>
{
  typedef Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 0.25))
          return {8.82e-02, 6.65e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 0.25 / 5.0))
          return {8.66e-03, 6.02e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {8.86e-02, 6.70e-02, 4.91e-02, 3.42e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1,
                                         Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1,
                                         Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

template class HyperbolicEocExpectations<Hyperbolic::ShockTubeTestCase<Dune::YaspGrid<1>, double>,
                                         Hyperbolic::ChooseDiscretizer::fv, 1,
                                         Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
