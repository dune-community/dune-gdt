// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_FINITEVOLUME_DEFAULT_HH
#define DUNE_GDT_PLAYGROUND_SPACES_FINITEVOLUME_DEFAULT_HH

#warning This header is deprecated, include <dune/gdt/playground/spaces/fv/default.hh> instead (19.11.2014)!
#include <dune/gdt/playground/spaces/fv/default.hh>

namespace Dune {
namespace GDT {
namespace Spaces {
namespace FiniteVolume {


template <class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DUNE_DEPRECATED_MSG("Use FV::Default instead (19.11.2014)!") Default
    : public FV::Default<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols>
{
public:
  template <class... Args>
  Default(Args&&... args)
    : FV::Default<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols>(std::forward<Args>(args)...)
  {
  }
};


} // namespace FiniteVolume
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_FINITEVOLUME_DEFAULT_HH
