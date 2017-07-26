// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH

#include <dune/xt/common/test/common.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/test/hyperbolic/eocstudy.hh>
#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems.hh"


template <class TestCaseType,
          Dune::GDT::NumericalFluxes numerical_flux,
          size_t reconstruction_order,
          Dune::GDT::TimeStepperMethods time_stepper,
          Dune::GDT::TimeStepperMethods rhs_time_stepper,
          Dune::GDT::TimeStepperSplittingMethods time_stepper_splitting>
struct hyperbolic_FV_discretization_base : public ::testing::Test
{
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
    DXTC_CONFIG.set("threading.partition_factor", 1u, true);
#if DXT_DISABLE_LARGE_TESTS
    TestCaseType test_case(/*num_refs = */ 1, /*divide_t_end_by = */ 5.0);
#else
    TestCaseType test_case;
#endif
    test_case.print_header(DXTC_LOG_INFO);
    DXTC_LOG_INFO << std::endl;
    typedef typename Hyperbolic::FvDiscretizer<TestCaseType,
                                               typename TestCaseType::GridType,
                                               double,
                                               TestCaseType::dimRange,
                                               TestCaseType::dimRangeCols,
                                               numerical_flux,
                                               reconstruction_order,
                                               time_stepper,
                                               rhs_time_stepper,
                                               time_stepper_splitting>
        Discretizer;
    Dune::GDT::Test::HyperbolicEocStudy<TestCaseType, Discretizer> eoc_study(test_case, {});
    XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DXTC_LOG_INFO));
  } // ... eoc_study()
}; // hyperbolic_FV_discretization_base

template <class TestCaseType>
struct hyperbolic_FV_discretization_entropy
    : public hyperbolic_FV_discretization_base<TestCaseType,
                                               Dune::GDT::NumericalFluxes::laxfriedrichs,
                                               0,
                                               Dune::GDT::TimeStepperMethods::explicit_rungekutta_second_order_ssp,
                                               Dune::GDT::TimeStepperMethods::matrix_exponential,
                                               Dune::GDT::TimeStepperSplittingMethods::strang>
{
}; // hyperbolic_FV_discretization_entropy

// template <class TestCaseType>
// struct hyperbolic_FV_discretization_godunovwithreconstruction_euler
//    : public hyperbolic_FV_discretization_base<TestCaseType,
//                                               Dune::GDT::NumericalFluxes::godunov_with_reconstruction,
//                                               Dune::GDT::TimeStepperMethods::explicit_euler>
//{
//}; // hyperbolic_FV_discretization_godunovwithreconstruction_euler

template <class TestCaseType>
struct hyperbolic_FV_discretization_laxfriedrichs_euler
    : public hyperbolic_FV_discretization_base<TestCaseType,
                                               Dune::GDT::NumericalFluxes::laxfriedrichs,
                                               0,
                                               Dune::GDT::TimeStepperMethods::explicit_euler,
                                               Dune::GDT::TimeStepperMethods::explicit_euler,
                                               Dune::GDT::TimeStepperSplittingMethods::fractional_step>
{
}; // hyperbolic_FV_discretization_godunov_laxfriedrichs_euler

template <class TestCaseType>
struct hyperbolic_FV_discretization_godunov_adaptiveRK
    : public hyperbolic_FV_discretization_base<TestCaseType,
                                               Dune::GDT::NumericalFluxes::godunov,
                                               0,
                                               Dune::GDT::TimeStepperMethods::dormand_prince,
                                               Dune::GDT::TimeStepperMethods::dormand_prince,
                                               Dune::GDT::TimeStepperSplittingMethods::fractional_step>
{
}; // hyperbolic_FV_discretization_godunov_adaptiveRK

#endif // DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
