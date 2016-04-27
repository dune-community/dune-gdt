// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#include "config.h"

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/grid.hh>

#include "../problems/spe10.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {9.48e-03, 2.64e-03};
    else if (type == "H1_semi")
      return {1.09e-01, 5.36e-02};
    else if (type == "energy")
      return {1.37e-01, 6.60e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.74e-03, 8.75e-04};
    else if (type == "H1_semi")
      return {4.91e-02, 2.26e-02};
    else if (type == "energy")
      return {5.60e-02, 2.58e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>, double,
                                                                        1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {9.43e-03, 2.55e-03};
    else if (type == "H1_semi")
      return {1.43e-01, 7.04e-02};
    else if (type == "energy")
      return {1.74e-01, 8.53e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>, double,
                                                                        1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.13e-03, 6.56e-04};
    else if (type == "H1_semi")
      return {6.41e-02, 3.24e-02};
    else if (type == "energy")
      return {7.56e-02, 3.81e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                 double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, conforming>,
                                                                                 double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;
template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                                 double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                                 double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID
