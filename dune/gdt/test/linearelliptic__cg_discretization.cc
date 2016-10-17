// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif

#include <dune/xt/common/test/main.hxx> // <- This one has to come first (includes the config.h)!

#include <dune/xt/common/test/common.hh>
#include <dune/xt/functions/spe10/model1.hh>
#include <dune/xt/la/container.hh>

#include <dune/gdt/spaces/interface.hh>

#include "linearelliptic/eocstudy.hh"
#include "linearelliptic/discretizers/cg.hh"

#include <dune/gdt/test/linearelliptic/problems.hh>
#include <dune/gdt/test/grids.hh>

struct linearelliptic_CG_discretization : public ::testing::Test
{
  typedef TESTCASETYPE TestCaseType;

  static void eoc_study()
  {
    using namespace Dune;
    using namespace Dune::GDT;
#if THIS_IS_A_BUILDBOT_BUILD
    TestCaseType test_case(/*num_refs=*/1); // As in: only 1!
#else
    TestCaseType test_case;
#endif
    test_case.print_header(DXTC_LOG_INFO);
    DXTC_LOG_INFO << std::endl;
    typedef LinearElliptic::CGDiscretizer<typename TestCaseType::GridType,
                                          XT::Grid::Layers::level,
                                          ChooseSpaceBackend::SPACE_BACKEND,
                                          XT::LA::Backends::LA_BACKEND,
                                          1,
                                          typename TestCaseType::ProblemType::RangeFieldType,
                                          1>
        Discretizer;
    Dune::GDT::Test::LinearEllipticEocStudy<TestCaseType, Discretizer> eoc_study(test_case);
    try {
      Dune::XT::Test::check_eoc_study_for_success(eoc_study, eoc_study.run(DXTC_LOG_INFO));
    } catch (Dune::XT::Common::Exceptions::spe10_data_file_missing&) {
      Dune::XT::Common::TimedLogger().get("gdt.test.linearelliptic.cg.discretization").warn()
          << "missing SPE10 data file!" << std::endl;
    }
  } // ... eoc_study()
}; // linearelliptic_CG_discretization

TEST_F(linearelliptic_CG_discretization, eoc_study)
{
  this->eoc_study();
}


namespace Dune {
namespace GDT {
namespace Test {

// YaspGrid<2, Dune::EquidistantOffsetCoordinates<double, 2>>
extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<Yasp2Grid, double, 1>, CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<Yasp2Grid, double, 1>, CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<Yasp2Grid, double, 1>, CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<Yasp2Grid, double, 1>, CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<Yasp2Grid, double, 1>, CG, 1>;


#if HAVE_ALUGRID


// ALUGrid< 2, 2, simplex, conforming >

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluConform2dGridType, double, 1>, CG,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluConform2dGridType, double, 1>, CG,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluConform2dGridType, double, 1>,
                                                    CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluConform2dGridType, double,
                                                                                          1>,
                                                    CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluConform2dGridType, double,
                                                                                        1>,
                                                    CG, 1>;


// ALUGrid< 2, 2, simplex, nonconforming >

extern template class LinearEllipticEocExpectations<LinearElliptic::AO2013TestCase<AluSimplex2dGridType, double, 1>, CG,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ER2007TestCase<AluSimplex2dGridType, double, 1>, CG,
                                                    1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::ESV2007TestCase<AluSimplex2dGridType, double, 1>,
                                                    CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::MixedBoundaryTestCase<AluSimplex2dGridType, double,
                                                                                          1>,
                                                    CG, 1>;

extern template class LinearEllipticEocExpectations<LinearElliptic::Spe10Model1TestCase<AluSimplex2dGridType, double,
                                                                                        1>,
                                                    CG, 1>;


#endif // HAVE_ALUGRID

} // namespace Test
} // namespace GDT
} // namespace Dune
