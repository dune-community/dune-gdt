// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

/**
  * This file is intended as a starting point for quick testing.
  */

#ifndef DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS
#define DUNE_STUFF_TEST_MAIN_CATCH_EXCEPTIONS 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_INFO_LOGGING 1
#endif
#ifndef DUNE_STUFF_TEST_MAIN_ENABLE_DEBUG_LOGGING
#define DUNE_STUFF_TEST_MAIN_ENABLE_DEBUG_LOGGING 1
#endif

#include <dune/stuff/test/main.hxx> // <- this one has to come first (includes the config.h)!

using namespace Dune;

TEST(empty, main)
{
}
