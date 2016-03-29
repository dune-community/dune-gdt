// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include "config.h"

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "problems/ER2007.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.11e-01, 6.09e-02, 1.65e-02};
    else if (type == "H1_semi" || type == "energy")
      return {3.66e-01, 2.98e-01, 1.46e-01};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {6.35e-02, 6.42e-03, 8.23e-04};
    else if (type == "H1_semi" || type == "energy")
      return {2.32e-01, 5.40e-02, 1.41e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {2.00e-01, 1.11e-01, 3.29e-02};
    else if (type == "H1_semi" || type == "energy")
      return {4.41e-01, 4.05e-01, 2.04e-01};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2, anything>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type)
  {
    if (type == "L2")
      return {1.51e-01, 1.27e-02, 1.50e-03};
    else if (type == "H1_semi" || type == "energy")
      return {3.69e-01, 8.28e-02, 2.12e-02};
    else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                            1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                            1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;
template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double, 1>,
                                             LinearElliptic::ChooseDiscretizer::swipdg, 2>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
