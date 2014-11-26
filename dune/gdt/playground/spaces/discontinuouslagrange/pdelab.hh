// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_PDELAB_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_PDELAB_HH

#warning This header is deprecated, include <dune/gdt/playground/spaces/dg/pdelab.hh> instead (21.11.2014)!
#include <dune/gdt/playground/spaces/dg/pdelab.hh>

namespace DiscontinuousLagrange {


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class
  DUNE_DEPRECATED_MSG("Use DG::PdelabBased instead (21.11.2014)!")
      PdelabBased
  : public DG::PdelabBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
public:
  template< class... Args >
  PdelabBased(Args&& ...args)
    : DG::PdelabBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >(std::forward< Args >(args)...)
  {}
};


} // namespace DiscontinuousLagrange

#endif // DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_PDELAB_HH
