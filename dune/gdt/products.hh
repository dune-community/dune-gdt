// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_PRODUCTS_HH
#define DUNE_GDT_PRODUCTS_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/bilinear-form.hh>


namespace Dune {
namespace GDT {


template <size_t r, // <- needs to be specified manually
          class GridViewType>
auto L2Product(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> weight = 1.,
               const int over_integrate = 0)
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto product = make_bilinear_form<r, 1, r, 1>(grid_view);
  product += LocalElementIntegralBilinearForm<E, r>(LocalProductIntegrand<E, r>(weight), over_integrate);
  return product;
}

template <class GridViewType>
auto L2Product(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
               const int over_integrate = 0)
{
  return L2Product<1>(grid_view, weight, over_integrate);
}


template <size_t r, // <- needs to be specified manually
          class GridViewType>
auto H1SemiProduct(const GridViewType& grid_view,
                   XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>,
                                               GridViewType::dimension,
                                               GridViewType::dimension> weight = 1.,
                   const int over_integrate = 0)
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto product = make_bilinear_form<r, 1, r, 1>(grid_view);
  product += LocalElementIntegralBilinearForm<E, r>(LocalLaplaceIntegrand<E, r>(weight), over_integrate);
  return product;
}

template <class GridViewType>
auto H1SemiProduct(const GridViewType& grid_view,
                   XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
                   const int over_integrate = 0)
{
  return H1SemiProduct<1>(grid_view, weight, over_integrate);
}


template <size_t r, // <- needs to be specified manually
          class GridViewType>
auto H1Product(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> l2_weight,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>,
                                           GridViewType::dimension,
                                           GridViewType::dimension> h1_semi_weight,
               const int over_integrate = 0)
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto product = make_bilinear_form<r, 1, r, 1>(grid_view);
  product += LocalElementIntegralBilinearForm<E, r>(LocalProductIntegrand<E, r>(l2_weight), over_integrate);
  product += LocalElementIntegralBilinearForm<E, r>(LocalLaplaceIntegrand<E, r>(h1_semi_weight), over_integrate);
  return product;
}

template <size_t r, // <- needs to be specified manually
          class GridViewType>
std::enable_if_t<r == GridViewType::dimension, BilinearForm<GridViewType, r, 1, r, 1>>
H1Product(const GridViewType& grid_view,
          XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>, r, r> weight = 1.,
          const int over_integrate = 0)
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto product = make_bilinear_form<r, 1, r, 1>(grid_view);
  product += LocalElementIntegralBilinearForm<E, r>(LocalProductIntegrand<E, r>(weight), over_integrate);
  product += LocalElementIntegralBilinearForm<E, r>(LocalLaplaceIntegrand<E, r>(weight), over_integrate);
  return product;
}

template <size_t r, // <- needs to be specified manually
          class GridViewType>
std::enable_if_t<r != GridViewType::dimension, BilinearForm<GridViewType, r, 1, r, 1>>
H1Product(const GridViewType& grid_view,
          XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
          const int over_integrate = 0)
{
  static_assert(XT::Grid::is_view<GridViewType>::value, "");
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto product = make_bilinear_form<r, 1, r, 1>(grid_view);
  product += LocalElementIntegralBilinearForm<E, r>(LocalProductIntegrand<E, r>(weight), over_integrate);
  product += LocalElementIntegralBilinearForm<E, r>(LocalLaplaceIntegrand<E, r>(weight), over_integrate);
  return product;
}

template <class GridViewType>
auto H1Product(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> l2_weight,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> h1_semi_weight,
               const int over_integrate = 0)
{
  return H1Product<1>(grid_view, l2_weight, h1_semi_weight, over_integrate);
}

template <class GridViewType>
auto H1Product(const GridViewType& grid_view,
               XT::Functions::GridFunction<XT::Grid::extract_entity_t<GridViewType>> weight = 1.,
               const int over_integrate = 0)
{
  return H1Product<1>(grid_view, weight, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_HH
