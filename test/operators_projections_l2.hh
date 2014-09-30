// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH

#include <dune/gdt/operators/projections.hh>

#include "operators_projections.hh"


template <class SpaceType>
struct L2ProjectionOperator
    : public ProjectionOperatorBase<SpaceType, Operators::L2Projection<typename SpaceType::GridViewType>>
{
};


#endif // DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
