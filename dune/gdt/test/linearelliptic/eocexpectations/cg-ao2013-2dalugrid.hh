// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)


#if HAVE_ALUGRID

#include <dune/grid/alugrid.hh>

#include "../problems/AO2013.hh"

#include <dune/gdt/test/linearelliptic/eocexpectations_base.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations

template <>
class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>,
                                    LinearElliptic::ChooseDiscretizer::cg, 1>
    : public internal::LinearEllipticEocExpectationsBase<1>
{
  typedef LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type); // ... results(...)
}; // LinearEllipticEocExpectations


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // HAVE_ALUGRID
