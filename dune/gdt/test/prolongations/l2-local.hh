// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH
#define DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH

#include <dune/gdt/prolongations/l2-local.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Test {


template <class SpaceType>
class L2LocalProlongationLocalizableOperatorTest
    : public LocalizableProlongationOperatorBase<SpaceType, SpaceType, L2LocalProlongationLocalizableOperator>
{
};


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_PROLONGATIONS_L2_LOCAL_HH
