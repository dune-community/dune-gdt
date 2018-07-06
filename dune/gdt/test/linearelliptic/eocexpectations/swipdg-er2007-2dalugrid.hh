// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_ER2007_2DALUGRID_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_ER2007_2DALUGRID_HH

#if HAVE_DUNE_ALUGRID

#include <dune/alugrid/grid.hh>

#include "../problems/ER2007.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {

// polorder 1, conforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, conforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    2> : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 1, nonconforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations

// polorder 2, nonconforming

template <>
class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    2> : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_ER2007_2DALUGRID_HH
