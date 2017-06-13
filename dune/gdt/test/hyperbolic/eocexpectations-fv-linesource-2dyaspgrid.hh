// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_POINTSOURCE_3DYASPGRID_HH
#define DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_POINTSOURCE_3DYASPGRID_HH

#include "config.h"

#include <dune/grid/yaspgrid.hh>

#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems/fokkerplanck/linesource.hh"
#include "eocexpectations.hh"


namespace Dune {
namespace GDT {
namespace Test {


template <>
class HyperbolicEocExpectations<Hyperbolic::LineSourceTestCase<Yasp2, double>,
                                Hyperbolic::ChooseDiscretizer::fv,
                                2,
                                NumericalFluxes::laxfriedrichs,
                                TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                TimeStepperMethods::matrix_exponential>
    : public internal::HyperbolicEocExpectationsBase<2>
{
  typedef Hyperbolic::LineSourceTestCase<Yasp2, double> TestCaseType;

public:
  static std::vector<double> results(const TestCaseType& test_case, const std::string type); // ... results(...)
}; // HyperbolicEocExpectations


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_HYPERBOLIC_EOCEXPECTATIONS_FV_POINTSOURCE_3DYASPGRID_HH
