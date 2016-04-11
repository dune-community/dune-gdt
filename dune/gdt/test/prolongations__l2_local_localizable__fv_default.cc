// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#include <dune/stuff/test/main.hxx>

#include "prolongations/l2-local.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types<SPACE_FV_SGRID_LEVEL(1, 1), SPACE_FV_SGRID_LEVEL(2, 1), SPACE_FV_SGRID_LEVEL(3, 1),
                       SPACE_FV_YASPGRID_LEVEL(1, 1), SPACE_FV_YASPGRID_LEVEL(2, 1), SPACE_FV_YASPGRID_LEVEL(3, 1)
#if HAVE_ALUGRID
                                                                                         ,
                       SPACE_FV_ALUCONFORMGRID_LEVEL(2, 1), SPACE_FV_ALUCONFORMGRID_LEVEL(3, 1),
                       SPACE_FV_ALUCUBEGRID_LEVEL(2, 1), SPACE_FV_ALUCUBEGRID_LEVEL(3, 1)
#endif // HAVE_ALUGRID
                       > SpaceTypes;

TYPED_TEST_CASE(L2LocalProlongationLocalizableOperatorTest, SpaceTypes);
TYPED_TEST(L2LocalProlongationLocalizableOperatorTest, constructible_by_ctor)
{
  this->constructible_by_ctor(1.45e-1);
}
TYPED_TEST(L2LocalProlongationLocalizableOperatorTest, constructible_by_factory)
{
  this->constructible_by_factory(1.45e-1);
}
TYPED_TEST(L2LocalProlongationLocalizableOperatorTest, produces_correct_results)
{
  this->produces_correct_results(1.45e-1);
}
