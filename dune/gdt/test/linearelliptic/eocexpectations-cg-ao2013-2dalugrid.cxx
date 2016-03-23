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


template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {1.17e-02, 5.32e-03};
      else
        return {1.30e-02, 6.81e-03, 1.96e-03, 4.47e-04};
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {3.63e-01, 2.97e-01};
      else
        return {4.11e-01, 3.53e-01, 1.90e-01, 8.78e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {9.22e-02, 5.04e-02};
      else
        return {9.76e-02, 5.97e-02, 3.20e-02, 1.54e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations

template <bool anything>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1, anything>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type)
  {
    if (type == "L2") {
      if (test_case.num_refinements() == 1)
        return {1.53e-02, 3.91e-03};
      else
        return {
            1.65e-02, 5.12e-03, 1.35e-03, 2.94e-04,
        };
    } else if (type == "H1_semi") {
      if (test_case.num_refinements() == 1)
        return {5.51e-01, 2.77e-01};
      else
        return {5.75e-01, 3.22e-01, 1.64e-01, 7.42e-02};
    } else if (type == "energy") {
      if (test_case.num_refinements() == 1)
        return {8.68e-02, 4.24e-02};
      else
        return {9.06e-02, 4.97e-02, 2.60e-02, 1.24e-02};
    } else
      EXPECT_TRUE(false) << "test results missing for type: " << type;
    return {};
  } // ... results(...)
}; // LinearEllipticEocExpectations


template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, conforming>, double,
                                                                            1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;
template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<ALUGrid<2, 2, simplex, nonconforming>,
                                                                            double, 1>,
                                             LinearElliptic::ChooseDiscretizer::cg, 1>;


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
