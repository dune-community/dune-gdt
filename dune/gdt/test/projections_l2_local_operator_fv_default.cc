// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "projections/l2-local.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_SGRID(1, 1), SPACE_FV_SGRID(2, 1), SPACE_FV_SGRID(3, 1), SPACE_FV_YASPGRID(1, 1),
                       SPACE_FV_YASPGRID(2, 1), SPACE_FV_YASPGRID(3, 1)
#if HAVE_ALUGRID && !defined(__GNUC__)
                                                    ,
                       SPACE_FV_ALUCONFORMGRID(2, 1), SPACE_FV_ALUCONFORMGRID(3, 1), SPACE_FV_ALUCUBEGRID(2, 1),
                       SPACE_FV_ALUCUBEGRID(3, 1)
#endif // HAVE_ALUGRID && !defined(__GNUC__)
                       > SpaceTypes;

TYPED_TEST_CASE(L2LocalProjectionOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProjectionOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor();
}
TYPED_TEST(L2LocalProjectionOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory();
}
TYPED_TEST(L2LocalProjectionOperatorTest, produces_correct_results)
{
  this->produces_correct_results(0.096226);
}
