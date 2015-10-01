// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PROLONGATIONS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_PROLONGATIONS_L2_HH

#include <dune/gdt/operators/prolongations/l2.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Tests {


template <class SpaceType>
class L2LocalProlongationLocalizableOperatorTest
    : public LocalizableProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationLocalizableOperator>
{
};


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_PROLONGATIONS_L2_HH
