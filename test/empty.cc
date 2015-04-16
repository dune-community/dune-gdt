// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

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
