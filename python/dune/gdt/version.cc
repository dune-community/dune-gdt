// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2021)

#include "config.h"

#include <dune/pybindxi/pybind11.h>

// see https://stackoverflow.com/questions/240353/convert-a-preprocessor-token-to-a-string
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)


PYBIND11_MODULE(_version, m)
{
  m.attr("__version__") = pybind11::str(TOSTRING(DUNE_GDT_VERSION));
}


#undef TOSTRING
#undef STRINGIFY
