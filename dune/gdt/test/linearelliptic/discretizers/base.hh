// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH


namespace Dune {
namespace GDT {
namespace LinearElliptic {


enum class ChooseDiscretizer
{
  cg,
  sipdg,
  swipdg,
  swipdg_affine_factor,
  swipdg_affine_tensor,
  block_ipdg
}; // enum class ChooseDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH
