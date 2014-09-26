// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_H1_HH
#define DUNE_GDT_PRODUCTS_H1_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/products/base.hh>
#include <dune/gdt/products/l2_internal.hh>
#include <dune/gdt/products/l2.hh>

#include "../localoperator/codim0.hh"
#include "../localevaluation/elliptic.hh"

#include "interfaces.hh"


namespace Dune {
namespace GDT {
namespace Products {

template <class FunctionImp>
using H1Evaluation = LocalEvaluation::Elliptic<FunctionImp, void /* = DiffusionTensorImp*/>;

template <class GridViewImp, class RangeImp, class SourceImp, class DerivedImp,
          template <class> class LocalEvaluationType>
using H1SemiLocalizableTraits =
    internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp, DerivedImp, LocalEvaluationType>;

template <class GridViewImp, class RangeImp, class SourceImp>
using H1SemiLocalizable = LocalizableForward<GridViewImp, RangeImp, SourceImp, H1SemiLocalizableTraits, H1Evaluation>;

template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class DerivedImp,
          template <class> class LocalEvaluationTemplate>
using H1SemiAssemblableTraits = internal::L2AssemblableTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                                              DerivedImp, LocalEvaluationTemplate>;
/**
 * \todo actual doc
 * \note this cannot be an alias because of the self-injection to base
 **/
template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
using H1SemiAssemblable =
    AssemblableForward<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, H1SemiAssemblableTraits, H1Evaluation>;

template <class GridViewImp, class FieldImp = double>
using H1SemiGeneric                         = ProductForward<GridViewImp, FieldImp, H1SemiLocalizable>;

} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_H1_HH
