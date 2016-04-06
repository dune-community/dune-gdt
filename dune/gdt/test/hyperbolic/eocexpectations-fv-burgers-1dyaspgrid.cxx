// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/burgers.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Tests {


template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperCombinations::godunov_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {8.96e-02, 3.87e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0/5.0))
          return {1.33e-02, 6.34e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {1.13e-01, 6.17e-02, 2.88e-02, 1.24e-02, 4.82e-03};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {8.64e-02, 3.72e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0/5.0))
          return {1.33e-02, 6.34e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {1.15e-01, 6.48e-02, 3.31e-02, 1.52e-02, 5.68e-03};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template< bool anything >
class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                 Hyperbolic::ChooseDiscretizer::fv,
                                 1,
                                 Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler,
                                 anything >
  : public internal::HyperbolicEocExpectationsBase< 1 >
{
  typedef Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 > TestCaseType;
public:
  static std::vector< double > results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L1") {
      if (test_case.num_refinements() == 1)
        if (DSC::FloatCmp::eq(test_case.t_end(), 1.0))
          return {1.03e-01, 5.58e-02};
        else if (DSC::FloatCmp::eq(test_case.t_end(), 1.0/5.0))
          return {1.71e-02, 7.63e-03};
        else
          EXPECT_TRUE(false) << "test results missing for t_end = " << DSC::toString(test_case.t_end());
      else
        return {1.76e-01, 1.30e-01, 7.79e-02, 3.85e-02, 1.41e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // HyperbolicEocExpectations

template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperCombinations::godunov_euler>;

template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperCombinations::godunov_adaptiveRK>;

template class HyperbolicEocExpectations< Hyperbolic::BurgersTestCase< Dune::YaspGrid< 1 >, double, 1 >,
                                          Hyperbolic::ChooseDiscretizer::fv,
                                          1,
                                          Hyperbolic::FluxTimeStepperCombinations::laxfriedrichs_euler>;


} // namespace Tests
} // namespace GDT
} // namespace Dune
