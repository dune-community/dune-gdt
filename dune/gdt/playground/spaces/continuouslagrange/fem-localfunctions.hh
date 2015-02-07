// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_PLAYGROUND_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH

#warning This header is deprecated, include <dune/gdt/playground/spaces/cg/fem-localfunctions.hh> instead (21.11.2014)!
#include <dune/gdt/playground/spaces/cg/fem-localfunctions.hh>

namespace Dune {
namespace GDT {
namespace Spaces {
namespace ContinuousLagrange {


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class
  DUNE_DEPRECATED_MSG("Use CG::FemLocalfunctionsBased instead (21.11.2014)!")
      FemLocalfunctionsBased
  : public CG::FemLocalfunctionsBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
public:
  template< class... Args >
  FemLocalfunctionsBased(Args&& ...args)
    : CG::FemLocalfunctionsBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >(std::forward< Args >(args)...)
  {}
};


} // namespace ContinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_CONTINUOUSLAGRANGE_FEM_LOCALFUNCTIONS_HH
