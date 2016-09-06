// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)


#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/AO2013.hh"
#include "../eocexpectations.hh"

namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg, 2>
    : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
