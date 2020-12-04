// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_NORMS_HH
#define DUNE_GDT_NORMS_HH

#include <dune/gdt/products.hh>


namespace Dune {
namespace GDT {


template <class GridViewType, size_t r>
double l2_norm(const GridViewType& grid_view,
               const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r>& function,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> weight = 1.,
               const int over_integrate = 0)
{
  return L2Product<r>(grid_view, weight, over_integrate).norm(function);
}

template <class GridViewType>
double l2_norm(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> function,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
               const int over_integrate = 0)
{
  return L2Product(grid_view, weight, over_integrate).norm(function);
}


template <class GridViewType, size_t r>
double h1_semi_norm(const GridViewType& grid_view,
                    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r>& function,
                    XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>,
                                                GridViewType::dimension,
                                                GridViewType::dimension> weight = 1.,
                    const int over_integrate = 0)
{
  return H1SemiProduct<r>(grid_view, weight, over_integrate).norm(function);
}

template <class GridViewType>
double h1_semi_norm(const GridViewType& grid_view,
                    XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> function,
                    XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>,
                                                GridViewType::dimension,
                                                GridViewType::dimension> weight = 1.,
                    const int over_integrate = 0)
{
  return H1SemiProduct<1>(grid_view, weight, over_integrate).norm(function);
}


template <class GridViewType, size_t r>
std::enable_if_t<r == GridViewType::dimension, double>
h1_norm(const GridViewType& grid_view,
        const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r>& function,
        XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> weight = 1.,
        const int over_integrate = 0)
{
  return H1Product<r>(grid_view, weight, over_integrate).norm(function);
}

template <class GridViewType, size_t r>
std::enable_if_t<r != GridViewType::dimension, double>
h1_norm(const GridViewType& grid_view,
        const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r>& function,
        XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
        const int over_integrate = 0)
{
  return H1Product<r>(grid_view, weight, over_integrate).norm(function);
}

template <class GridViewType, size_t r>
double h1_norm(const GridViewType& grid_view,
               const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r>& function,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> l2_weight,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>,
                                           GridViewType::dimension,
                                           GridViewType::dimension> h1_semi_weight,
               const int over_integrate = 0)
{
  return H1Product<r>(grid_view, l2_weight, h1_semi_weight, over_integrate).norm(function);
}

template <class GridViewType>
double h1_norm(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> function,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
               const int over_integrate = 0)
{
  return H1Product(grid_view, weight, over_integrate).norm(function);
}

template <class GridViewType>
double h1_norm(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> function,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> l2_weight,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> h1_semi_weight,
               const int over_integrate = 0)
{
  return H1Product(grid_view, l2_weight, h1_semi_weight, over_integrate).norm(function);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_NORMS_HH
