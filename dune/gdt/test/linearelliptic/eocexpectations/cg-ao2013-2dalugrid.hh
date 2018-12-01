// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_AO2013_2DALUGRID_HH
#define DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_AO2013_2DALUGRID_HH

#if HAVE_DUNE_ALUGRID

#  include <dune/alugrid/grid.hh>

#  include "../problems/AO2013.hh"

#  include <dune/gdt/test/linearelliptic/eocexpectations_base.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations


template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluCube2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg,
                                    1> : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluCube2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type);
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_ALUGRID

#endif // DUNE_GDT_TEST_LINEARELLIPTIC_EOCEXPECTATIONS_CG_AO2013_2DALUGRID_HH
