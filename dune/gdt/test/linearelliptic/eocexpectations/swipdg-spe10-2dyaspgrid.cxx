// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include "../problems/spe10.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::
                                        Spe10Model1TestCase<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double,
                                                                                                                 2>>,
                                                            double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>, double,
                                              1>
      TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.73e-02, 2.35e-02};
    else if (type == "H1_semi")
      return {4.64e-01, 6.31e-01};
    else if (type == "energy")
      return {1.03e+00, 1.49e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::
                                        Spe10Model1TestCase<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double,
                                                                                                                 2>>,
                                                            double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::Spe10Model1TestCase<Dune::YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>, double,
                                              1>
      TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.82e-02, 7.50e-03};
    else if (type == "H1_semi")
      return {5.53e-01, 4.55e-01};
    else if (type == "energy")
      return {1.62e+00, 1.33e+00};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      Spe10Model1TestCase<Dune::YaspGrid<2,
                                                                         Dune::EquidistantOffsetCoordinates<double, 2>>,
                                                          double, 1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class
    LinearEllipticEocExpectations<LinearElliptic::
                                      Spe10Model1TestCase<Dune::YaspGrid<2,
                                                                         Dune::EquidistantOffsetCoordinates<double, 2>>,
                                                          double, 1>,
                                  LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune
