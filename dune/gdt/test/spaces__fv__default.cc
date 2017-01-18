// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as  BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014, 2016 - 2017)
//   Rene Milk       (2014, 2016)

// This one has to come first (includes the config.h)!
#include <dune/xt/common/test/main.hxx>

#include "spaces/fv.hh"
#include "spaces/fv/default.hh"


typedef testing::Types<SPACES_FV
#if HAVE_ALUGRID
                       ,
                       SPACES_FV_ALUGRID
#endif
                       >
    FV_Spaces;

TYPED_TEST_CASE(FV_Space, FV_Spaces);
TYPED_TEST(FV_Space, fulfills_interface)
{
  this->fulfills_interface();
}
TYPED_TEST(FV_Space, mapper_fulfills_interface)
{
  this->mapper_fulfills_interface();
}
TYPED_TEST(FV_Space, basefunctionset_fulfills_interface)
{
  this->basefunctionset_fulfills_interface();
}
TYPED_TEST(FV_Space, check_for_correct_copy)
{
  this->check_for_correct_copy();
}
