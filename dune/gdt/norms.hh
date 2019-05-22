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

#include <dune/grid/common/gridview.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/localizable-bilinear-form.hh>


namespace Dune {
namespace GDT {


template <class GridViewType, size_t r, class SF, class RF, class F = double>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value,
                 LocalizableBilinearFormBase<GridViewType, r, 1, SF, F, r, 1, RF>>
make_localizable_l2_product(
    const GridViewType& grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, SF>& left,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, RF>& right,
    const int over_integrate = 0)
{
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto localizable_product = make_localizable_bilinear_form(grid_view, left, right);
  localizable_product.append(std::make_unique<LocalElementIntegralBilinearForm<E, r, 1, RF, F, r, 1, SF>>(
      LocalElementProductIntegrand<E, r, RF, F, SF>(), over_integrate));
  return localizable_product;
}


template <class GridViewType, size_t r, class SF, class RF, class F = double>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value, F>
l2_product(const GridViewType& grid_view,
           const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, SF>& left,
           const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, RF>& right,
           const int over_integrate = 0)
{
  auto product = make_localizable_l2_product(grid_view, left, right, over_integrate);
  product.assemble();
  return product.result();
}


template <class GridViewType, size_t r, class R, class F = double>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value, F>
l2_norm(const GridViewType& grid_view,
        const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, R>& func,
        const int over_integrate = 0)
{
  return std::sqrt(l2_product(grid_view, func, func, over_integrate));
}


template <class GridViewType, class F, size_t r>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value, LocalizableBilinearFormBase<GridViewType, r, 1, F, F, 1, 1, F>>
make_localizable_elliptic_product(
    const GridViewType& grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>,
                                               GridViewType::dimension,
                                               GridViewType::dimension,
                                               F>& diffusion_tensor,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, F>& left,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, F>& right,
    const int over_integrate = 0)
{
  using E = XT::Grid::extract_entity_t<GridViewType>;
  auto localizable_product = make_localizable_bilinear_form(grid_view, left, right);
  localizable_product.append(LocalElementIntegralBilinearForm<E, r, 1, F, F, r, 1, F>(
      LocalLaplaceIntegrand<E, r, F>(diffusion_tensor), over_integrate));
  return localizable_product;
}


template <class GridViewType, class F, size_t r>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value, F>
elliptic_product(const GridViewType& grid_view,
                 const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>,
                                                            GridViewType::dimension,
                                                            GridViewType::dimension,
                                                            F>& diffusion_tensor,
                 const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, F>& left,
                 const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, F>& right,
                 const int over_integrate = 0)
{
  auto product = make_localizable_elliptic_product(grid_view, diffusion_tensor, left, right, over_integrate);
  product.assemble();
  return product.result();
}


template <class GridViewType, class F, size_t r>
std::enable_if_t<XT::Grid::is_view<GridViewType>::value, F>
elliptic_norm(const GridViewType& grid_view,
              const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>,
                                                         GridViewType::dimension,
                                                         GridViewType::dimension,
                                                         F>& diffusion_tensor,
              const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<GridViewType>, r, 1, F>& func,
              const int over_integrate = 0)
{
  return std::sqrt(elliptic_product(grid_view, diffusion_tensor, func, func, over_integrate));
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_NORMS_HH
