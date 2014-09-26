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
using H1Evaluation                          = LocalEvaluation::Elliptic<FunctionImp, void /* = DiffusionTensorImp*/>;
template <class GridViewImp, class FieldImp = double>
using H1SemiTraits                          = internal::L2BaseTraits<GridViewImp, FieldImp, H1Evaluation>;

template <class GridViewImp, class FieldImp>
using H1SemiBase = internal::L2Base<GridViewImp, FieldImp, H1Evaluation>;

template <class GridViewImp, class RangeImp, class SourceImp>
struct H1SemiLocalizable;

template <class GridViewImp, class RangeImp, class SourceImp, class DerivedImp,
          template <class> class LocalEvaluationType>
using H1SemiLocalizableTraits =
    internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp, DerivedImp, LocalEvaluationType>;

/**
 * \todo actual doc
 * \note this cannot be an alias because of the self-injection to base
 **/
template <class GridViewImp, class RangeImp, class SourceImp>
struct H1SemiLocalizable
    : public LocalizableForward<GridViewImp, RangeImp, SourceImp, H1SemiLocalizable<GridViewImp, RangeImp, SourceImp>,
                                H1SemiLocalizableTraits, H1Evaluation>
{
  typedef LocalizableForward<GridViewImp, RangeImp, SourceImp, H1SemiLocalizable<GridViewImp, RangeImp, SourceImp>,
                             H1SemiLocalizableTraits, H1Evaluation> BaseType;
  template <class... Args>
  explicit H1SemiLocalizable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};

template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class DerivedImp,
          template <class> class LocalEvaluationTemplate>
using H1SemiAssemblableTraits = internal::L2AssemblableTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                                              DerivedImp, LocalEvaluationTemplate>;
/**
 * \todo actual doc
 * \note this cannot be an alias because of the self-injection to base
 **/
template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
struct H1SemiAssemblable
    : public AssemblableForward<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                H1SemiAssemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>,
                                H1SemiAssemblableTraits, H1Evaluation>
{
  typedef AssemblableForward<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                             H1SemiAssemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>,
                             H1SemiAssemblableTraits, H1Evaluation> BaseType;
  template <class... Args>
  explicit H1SemiAssemblable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};

template <class GridViewImp, class FieldImp = double>
class H1SemiGeneric;


template <class GridViewImp, class FieldImp = double>
class H1SemiGenericTraits
{
public:
  typedef H1SemiGeneric<GridViewImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
};


template <class GridViewImp, class FieldImp>
class H1SemiGeneric : public ProductInterface<H1SemiGenericTraits<GridViewImp, FieldImp>>
{
public:
  typedef H1SemiGenericTraits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  H1SemiGeneric(const GridViewType& grd_vw, const size_t over_integrate = 0)
    : grid_view_(grd_vw)
    , over_integrate_(over_integrate)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType apply2(
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <int dimRangeRows, int dimRangeCols>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& range,
                   const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& source) const
  {
    typedef Stuff::
        LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols>
            FunctionType;
    H1SemiLocalizable<GridViewType, FunctionType, FunctionType> product_operator(
        grid_view_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class H1SemiGeneric


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_H1_HH
