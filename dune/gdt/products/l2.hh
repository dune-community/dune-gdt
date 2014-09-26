// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_HH
#define DUNE_GDT_PRODUCTS_L2_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/products/l2_internal.hh>
#include <dune/gdt/products/generic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/products/base.hh>

namespace Dune {
namespace GDT {
namespace Products {

template <class GridViewImp, class RangeImp, class SourceImp>
using L2Localizable =
    GenericLocalizable<GridViewImp, RangeImp, SourceImp, internal::L2LocalizableTraits, LocalEvaluation::Product>;

template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
using L2Assemblable = GenericAssemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                         internal::L2AssemblableTraits, LocalEvaluation::Product>;

template <class GridViewImp, class FieldImp = double>
using L2                                    = GenericProduct<GridViewImp, FieldImp, L2Localizable>;

} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_HH
