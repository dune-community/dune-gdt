// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

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
  swipdg_affine_tensor
}; // enum class ChooseDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH
