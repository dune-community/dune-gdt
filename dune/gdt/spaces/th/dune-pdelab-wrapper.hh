// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_TH_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_SPACES_TH_DUNE_PDELAB_WRAPPER_HH

#include <dune/gdt/spaces/product.hh>

namespace Dune {
namespace GDT {


template <class GridViewImp, int polynomialOrder, size_t domainDim, class RangeFieldImp>
class DunePdelabTaylorHoodSpaceWrapper
    : public DefaultProductSpace<DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1>,
                                 DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1>>
{
public:
  typedef DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1> VelocitySpaceType;
  typedef DunePdelabCgSpaceWrapper<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1> PressureSpaceType;
  typedef DefaultProductSpace<VelocitySpaceType, PressureSpaceType> BaseType;
  DunePdelabTaylorHoodSpaceWrapper(const GridViewImp& grid_view)
    : BaseType(VelocitySpaceType(grid_view), PressureSpaceType(grid_view))
  {
  }
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_TH_DUNE_PDELAB_WRAPPER_HH
