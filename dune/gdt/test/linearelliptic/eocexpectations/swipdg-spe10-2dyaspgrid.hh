// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_SPE10_2DYASPGRID_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_SPE10_2DYASPGRID_HH

#include <dune/grid/yaspgrid.hh>

#include "../problems/spe10.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations

template <>
class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::swipdg,
                                    2> : public internal::LinearEllipticEocExpectationsBase<2>
{
  typedef LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_SWIPDG_SPE10_2DYASPGRID_HH
