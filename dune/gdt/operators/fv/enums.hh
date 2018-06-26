// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_ENUMS_HH
#define DUNE_GDT_OPERATORS_FV_ENUMS_HH

namespace Dune {
namespace GDT {


enum class NumericalFluxes
{
  force,
  godunov,
  kinetic,
  laxfriedrichs,
  laxwendroff,
  local_laxfriedrichs,
  musta
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENUMS_HH
