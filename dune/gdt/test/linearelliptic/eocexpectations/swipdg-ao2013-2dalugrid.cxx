// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "problems/AO2013.hh"
#include "eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {5.33e-02, 1.69e-02};
      else
        return {5.57e-02, 1.99e-02, 5.54e-03, 1.29e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {3.82e-01, 2.29e-01};
      else
        return {4.32e-01, 2.93e-01, 1.50e-01, 6.54e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {2.32e-01, 1.15e-01};
      else
        return {2.53e-01, 1.46e-01, 7.45e-02, 3.35e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {1.18e-02, 2.12e-03};
      else
        return {1.18e-02, 2.11e-03, 3.89e-04, 7.76e-05};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {1.67e-01, 5.58e-02};
      else
        return {1.69e-01, 5.96e-02, 1.94e-02, 6.04e-03};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {6.96e-02, 2.59e-02};
      else
        return {7.41e-02, 3.36e-02, 1.62e-02, 7.03e-03};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {9.88e-02, 3.08e-02};
      else
        return {1.05e-01, 3.90e-02, 1.27e-02, 3.13e-03};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {5.95e-01, 3.31e-01};
      else
        return {6.52e-01, 4.24e-01, 2.20e-01, 9.59e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {3.39e-01, 1.63e-01};
      else
        return {3.67e-01, 2.08e-01, 1.06e-01, 4.75e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {2.77e-02, 5.49e-03};
      else
        return {2.77e-02, 5.53e-03, 8.36e-04, 1.29e-04};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {2.74e-01, 1.13e-01};
      else
        return {2.75e-01, 1.17e-01, 3.67e-02, 1.10e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {1.11e-01, 4.63e-02};
      else
        return {1.15e-01, 5.29e-02, 2.34e-02, 1.05e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                            1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                            1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;
template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
