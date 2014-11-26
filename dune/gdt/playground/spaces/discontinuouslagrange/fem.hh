// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH

#warning This header is deprecated, include <dune/gdt/playground/spaces/dg/fem.hh> instead (21.11.2014)!
#include <dune/gdt/playground/spaces/dg/fem.hh>

namespace Dune {
namespace GDT {
namespace Spaces {
namespace DiscontinuousLagrange {


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class
  DUNE_DEPRECATED_MSG("Use DG::FemBased instead (21.11.2014)!")
      FemBased
  : public DG::FemBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
public:
  template< class... Args >
  FemBased(Args&& ...args)
    : DG::FemBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >(std::forward< Args >(args)...)
  {}
};


} // namespace DiscontinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_DISCONTINUOUSLAGRANGE_FEM_HH
