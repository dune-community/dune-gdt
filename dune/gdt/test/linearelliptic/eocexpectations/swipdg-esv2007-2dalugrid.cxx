// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/ESV2007.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.82e-02, 4.53e-03, 1.12e-03, 2.78e-04};
    else if (type == "H1_semi" || type == "energy")
      return {1.48e-01, 7.28e-02, 3.62e-02, 1.80e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {8.55e-04, 1.06e-04, 1.31e-05, 1.63e-06};
    else if (type == "H1_semi" || type == "energy")
      return {1.41e-02, 3.56e-03, 8.91e-04, 2.23e-04};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.31e-02, 5.96e-03, 1.50e-03, 3.76e-04};
    else if (type == "H1_semi" || type == "energy")
      return {1.49e-01, 7.33e-02, 3.63e-02, 1.81e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {6.78e-04, 8.23e-05, 1.02e-05, 1.26e-06};
    else if (type == "H1_semi" || type == "energy")
      return {1.32e-02, 3.30e-03, 8.25e-04, 2.06e-04};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                             1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                             1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;
template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                             double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                             double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
