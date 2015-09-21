// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH

#include <dune/gdt/operators/projections/l2.hh>

#include "../projections.hh"

namespace Dune {
namespace GDT {
namespace Tests {


template< class SpaceType >
class L2LocalProjectionLocalizableOperatorTest
  : public LocalizableProjectionOperatorBase< SpaceType, L2LocalProjectionLocalizableOperator<
        typename SpaceType::GridViewType,
        typename internal::OperatorBaseTraits< SpaceType >::FunctionType,
        typename internal::OperatorBaseTraits< SpaceType >::DiscreteFunctionType > >
{};


template< class SpaceType >
class L2LocalProjectionOperatorTest
  : public ProjectionOperatorBase< SpaceType
                                 , L2LocalProjectionOperator< typename SpaceType::GridViewType, double > >
{};


} // namespace Tests
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_OPERATORS_PROJECTIONS_L2_HH
