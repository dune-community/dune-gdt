// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_PROJECTIONS_DIRICHLET_HH
#define DUNE_GDT_PROJECTIONS_DIRICHLET_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/intersection.hh>
#include <dune/xt/grid/layers.hh>

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
  typedef XT::Grid::BoundaryInfo<IntersectionType> BoundaryInfoType;

  template <class... Args>
  explicit DirichletProjectionLocalizableOperator(const BoundaryInfoType& boundary_info, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_operator_(boundary_info)
  {
    this->append(local_operator_);
    this->range().vector() *= 0.0;
  }

private:
  const LocalDirichletProjectionOperator<BoundaryInfoType> local_operator_;
}; // class DirichletProjectionLocalizableOperator


template <class GridViewType, class SourceType, class RangeType>
std::unique_ptr<DirichletProjectionLocalizableOperator<GridViewType, SourceType, RangeType>>
make_localizable_dirichlet_projection_operator(
    const GridViewType& grid_view,
    const XT::Grid::BoundaryInfo<typename XT::Grid::Intersection<GridViewType>::Type>& boundary_info,
    const SourceType& source, RangeType& range)
{
  return Dune::XT::Common::make_unique<DirichletProjectionLocalizableOperator<GridViewType, SourceType, RangeType>>(
      boundary_info, grid_view, source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_DIRICHLET_HH
