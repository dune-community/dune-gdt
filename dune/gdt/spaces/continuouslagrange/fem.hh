// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH

#warning This header is deprecated, include <dune/gdt/spaces/cg/fem.hh> instead (21.11.2014)!
#include <dune/gdt/spaces/cg/fem.hh>

namespace ContinuousLagrange {


template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class
  DUNE_DEPRECATED_MSG("Use CG::FemBased instead (21.11.2014)!")
      FemBased
  : public CG::FemBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
public:
  template< class... Args >
  FemBased(Args&& ...args)
    : CG::FemBased< GridPartImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >(std::forward< Args >(args)...)
  {}
};


} // namespace ContinuousLagrange

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_FEM_HH
