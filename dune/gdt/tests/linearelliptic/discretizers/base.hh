// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH
#define DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH


namespace Dune {
namespace GDT {
namespace LinearElliptic {


enum class ChooseDiscretizer
{
    cg
  , sipdg
  , swipdg
}; // enum class ChooseDiscretizer


} // namespace LinearElliptic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TESTS_LINEARELLIPTIC_DISCRETIZERS_BASE_HH
