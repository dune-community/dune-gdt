// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_ESV2007_2DALUGRID_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_ESV2007_2DALUGRID_HH

#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/ESV2007.hh"
#include "../eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type);
}; // LinearEllipticEocExpectations

template <>
class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& /*test_case*/, const std::string type);
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_ESV2007_2DYASPGRID_HH
