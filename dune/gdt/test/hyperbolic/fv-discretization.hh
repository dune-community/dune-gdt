// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
#define DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH

#ifndef THIS_IS_A_BUILDBOT_BUILD
#define THIS_IS_A_BUILDBOT_BUILD 0
#endif

#include <dune/stuff/test/common.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/test/hyperbolic/eocstudy.hh>
#include <dune/gdt/test/hyperbolic/discretizers/fv.hh>

#include "problems.hh"


template <class TestCaseType, Dune::GDT::NumericalFluxes numerical_flux, Dune::GDT::TimeStepperMethods time_stepper>
struct hyperbolic_FV_discretization_base : public ::testing::Test
{
  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
#if THIS_IS_A_BUILDBOT_BUILD
    TestCaseType test_case(/*num_refs = */ 1, /*divide_t_end_by = */ 5.0);
#else
    TestCaseType test_case;
#endif
    test_case.print_header(DSC_LOG_INFO);
    DSC_LOG_INFO << std::endl;
    typedef typename Hyperbolic::FvDiscretizer<TestCaseType,
                                               typename TestCaseType::GridType,
                                               double,
                                               TestCaseType::dimRange,
                                               TestCaseType::dimRangeCols,
                                               numerical_flux,
                                               time_stepper>
        Discretizer;
    Tests::HyperbolicEocStudy<TestCaseType, Discretizer> eoc_study(test_case, {});
    Stuff::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DSC_LOG_INFO));
  } // ... eoc_study()
}; // hyperbolic_FV_discretization_base

template <class TestCaseType>
struct hyperbolic_FV_discretization_godunov_euler
    : public hyperbolic_FV_discretization_base<TestCaseType, Dune::GDT::NumericalFluxes::godunov,
                                               Dune::GDT::TimeStepperMethods::explicit_euler>
{
}; // hyperbolic_FV_discretization_godunov_euler

template <class TestCaseType>
struct hyperbolic_FV_discretization_godunovwithreconstruction_euler
    : public hyperbolic_FV_discretization_base<TestCaseType, Dune::GDT::NumericalFluxes::godunov_with_reconstruction,
                                               Dune::GDT::TimeStepperMethods::explicit_euler>
{
}; // hyperbolic_FV_discretization_godunovwithreconstruction_euler

template <class TestCaseType>
struct hyperbolic_FV_discretization_laxfriedrichs_euler
    : public hyperbolic_FV_discretization_base<TestCaseType, Dune::GDT::NumericalFluxes::laxfriedrichs,
                                               Dune::GDT::TimeStepperMethods::explicit_euler>
{
}; // hyperbolic_FV_discretization_godunov_laxfriedrichs_euler

template <class TestCaseType>
struct hyperbolic_FV_discretization_godunov_adaptiveRK
    : public hyperbolic_FV_discretization_base<TestCaseType, Dune::GDT::NumericalFluxes::godunov,
                                               Dune::GDT::TimeStepperMethods::dormand_prince>
{
}; // hyperbolic_FV_discretization_godunov_adaptiveRK

#endif // DUNE_GDT_TEST_HYPERBOLIC_FV_DISCRETIZATION_HH
