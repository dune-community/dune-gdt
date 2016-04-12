// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PROJECTIONS_DIRICHLET_HH
#define DUNE_GDT_PROJECTIONS_DIRICHLET_HH

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/layers.hh>

#include <dune/gdt/local/operators/dirichlet-projection.hh>
#include <dune/gdt/operators/base.hh>

namespace Dune {
namespace GDT {


template <class GridViewImp, class SourceImp, class RangeImp, class FieldImp = double>
class DirichletProjectionLocalizableOperator : public LocalizableOperatorBase<GridViewImp, SourceImp, RangeImp>
{
  typedef LocalizableOperatorBase<GridViewImp, SourceImp, RangeImp> BaseType;

public:
  using typename BaseType::IntersectionType;
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionType> BoundaryInfoType;

  template <class... Args>
  explicit DirichletProjectionLocalizableOperator(const BoundaryInfoType& boundary_info, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_operator_(boundary_info)
  {
    this->add(local_operator_);
    this->range().vector() *= 0.0;
  }

private:
  const LocalDirichletProjectionOperator<BoundaryInfoType> local_operator_;
}; // class DirichletProjectionLocalizableOperator


template <class GridViewType, class SourceType, class RangeType>
std::unique_ptr<DirichletProjectionLocalizableOperator<GridViewType, SourceType, RangeType>>
make_localizable_dirichlet_projection_operator(
    const GridViewType& grid_view,
    const Stuff::Grid::BoundaryInfoInterface<typename Stuff::Grid::Intersection<GridViewType>::Type>& boundary_info,
    const SourceType& source, RangeType& range)
{
  return DSC::make_unique<DirichletProjectionLocalizableOperator<GridViewType, SourceType, RangeType>>(
      boundary_info, grid_view, source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_DIRICHLET_HH
