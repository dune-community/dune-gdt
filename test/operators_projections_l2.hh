// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH

#include <dune/gdt/operators/projections.hh>

#include "operators_projections.hh"


template< class SpaceType >
class L2ProjectionOperator
  : public ProjectionOperatorBase< SpaceType, Dune::GDT::Operators::L2Projection< typename SpaceType::GridViewType > >
{
  typedef ProjectionOperatorBase< SpaceType, Dune::GDT::Operators::L2Projection< typename SpaceType::GridViewType > >
      BaseType;
public:
  using typename BaseType::RangeFieldType;

  void free_project_l2_function_works(const RangeFieldType& tolerance = 1e-15)
  {
    this->vector_ *= 0.0;
    Dune::GDT::project_l2(this->function_, this->discrete_function_);
    this->measure_error(tolerance);
  }
};


#endif // DUNE_GDT_OPERATORS_PROJECTIONS_L2_HH
