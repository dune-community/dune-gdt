// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1

#include <dune/stuff/test/main.hxx>

#include "prolongations/l2.hh"
#include "spaces/fv/default.hh"

using namespace Dune::GDT::Test;


typedef testing::Types< SPACE_FV_SGRID_LEVEL(1, 1)
                      , SPACE_FV_SGRID_LEVEL(2, 1)
                      , SPACE_FV_SGRID_LEVEL(3, 1)
                      , SPACE_FV_YASPGRID_LEVEL(1, 1)
                      , SPACE_FV_YASPGRID_LEVEL(2, 1)
                      , SPACE_FV_YASPGRID_LEVEL(3, 1)
#if HAVE_ALUGRID && !defined(__GNUC__)
                      , SPACE_FV_ALUCONFORMGRID_LEVEL(2, 1)
                      , SPACE_FV_ALUCONFORMGRID_LEVEL(3, 1)
                      , SPACE_FV_ALUCUBEGRID_LEVEL(2, 1)
                      , SPACE_FV_ALUCUBEGRID_LEVEL(3, 1)
#endif // HAVE_ALUGRID && !defined(__GNUC__)
                      > SpaceTypes;

TYPED_TEST_CASE(L2ProlongationOperatorTest, SpaceTypes);
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_ctor) {
  this->constructible_by_ctor(1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, constructible_by_factory) {
  this->constructible_by_factory(1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, produces_correct_results) {
  this->produces_correct_results(1.45e-1);
}
TYPED_TEST(L2ProlongationOperatorTest, free_function_callable) {
  this->free_function_callable(1.45e-1);
}
