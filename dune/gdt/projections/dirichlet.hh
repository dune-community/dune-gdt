// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_PROJECTIONS_DIRICHLET_HH
#define DUNE_GDT_PROJECTIONS_DIRICHLET_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/layers.hh>

#include <dune/gdt/local/operators/dirichlet-projection.hh>
#include <dune/gdt/operators/base.hh>

namespace Dune {
namespace GDT {


template <class GridLayerImp, class SourceImp, class RangeImp, class FieldImp = double>
class DirichletProjectionLocalizableOperator : public LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp>
{
  typedef LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp> BaseType;

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


template <class GridLayerType, class SourceType, class RangeType>
std::unique_ptr<DirichletProjectionLocalizableOperator<GridLayerType, SourceType, RangeType>>
make_localizable_dirichlet_projection_operator(
    const GridLayerType& grid_layer,
    const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info,
    const SourceType& source,
    RangeType& range)
{
  return Dune::XT::Common::make_unique<DirichletProjectionLocalizableOperator<GridLayerType, SourceType, RangeType>>(
      boundary_info, grid_layer, source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_DIRICHLET_HH
