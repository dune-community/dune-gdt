// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_HH
#define DUNE_GDT_PRODUCTS_L2_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/grid/common/gridview.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/products/l2_internal.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/products/base.hh>

namespace Dune {
namespace GDT {
namespace Products {

template< class GridViewImp, class RangeImp, class SourceImp, class AliasedType,
          template <class, class, class, class, template <class> class> class TraitsTemplate,
          template <class> class LocalEvaluationTemplate>
class LocalizableForward
  : public LocalizableBase< TraitsTemplate< GridViewImp, RangeImp, SourceImp, AliasedType, LocalEvaluationTemplate > >
  , public internal::L2Base< GridViewImp, typename RangeImp::RangeFieldType, LocalEvaluationTemplate >
{
public:
  typedef TraitsTemplate< GridViewImp, RangeImp, SourceImp, AliasedType, LocalEvaluationTemplate > Traits;
private:
  typedef Products::LocalizableBase< Traits > LocalizableBaseType;
  typedef internal::L2Base< GridViewImp, typename RangeImp::RangeFieldType, LocalEvaluationTemplate> L2BaseType;
public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  LocalizableForward(const GridViewType& grd_vw, const RangeType& rng, const SourceType& src,
                const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, src)
    , L2BaseType(over_integrate)
  {}

  LocalizableForward(const GridViewType& grd_vw, const RangeType& rng, const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, rng)
    , L2BaseType(over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Localizable

/**
 * \todo actual doc
 * \note this cannot be an alias because of the self-injection to base
 **/
template< class GridViewImp, class RangeImp, class SourceImp >
struct L2Localizable
    : public LocalizableForward<GridViewImp, RangeImp, SourceImp, L2Localizable<GridViewImp, RangeImp, SourceImp>,
                           internal::L2LocalizableTraits, LocalEvaluation::Product> {
  typedef LocalizableForward<GridViewImp, RangeImp, SourceImp, L2Localizable<GridViewImp, RangeImp, SourceImp>,
                        internal::L2LocalizableTraits, LocalEvaluation::Product> BaseType;
  template <class... Args>
  explicit L2Localizable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};

template< class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class AliasedType,
          template <class, class, class, class, class, template <class> class > class TraitsTemplate,
          template <class> class LocalEvaluationTemplate>
class AssemblableForward
  : public AssemblableBase< TraitsTemplate< MatrixImp, RangeSpaceImp, GridViewImp,
                                            SourceSpaceImp, AliasedType, LocalEvaluationTemplate > >
  , public internal::L2Base< GridViewImp, typename RangeSpaceImp::RangeFieldType, LocalEvaluationTemplate >
{
public:
  typedef TraitsTemplate< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                          AliasedType, LocalEvaluationTemplate> Traits;
private:
  typedef Products::AssemblableBase< Traits >
    AssemblableBaseType;
  typedef internal::L2Base< GridViewImp, typename RangeSpaceImp::RangeFieldType, LocalEvaluationTemplate> L2BaseType;
public:
  typedef typename Traits::GridViewType     GridViewType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::MatrixType       MatrixType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  using AssemblableBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  AssemblableForward(MatrixType& mtrx,
                const RangeSpaceType& rng_scp,
                const GridViewType& grd_vw,
                const SourceSpaceType& src_scp,
                const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_scp, grd_vw, src_scp)
    , L2BaseType(over_integrate)
  {}

  AssemblableForward(MatrixType& mtrx,
                const RangeSpaceType& rng_scp,
                const GridViewType& grd_vw,
                const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_scp, grd_vw, rng_scp)
    , L2BaseType(over_integrate)
  {}

  AssemblableForward(MatrixType& matrix,
                const RangeSpaceType& range_space,
                const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , L2BaseType(over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class AssemblableForward

/**
 * \todo actual doc
 * \note this cannot be an alias because of the self-injection to base
 **/
template< class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
struct L2Assemblable
    : public AssemblableForward<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                                L2Assemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>,
                           internal::L2AssemblableTraits, LocalEvaluation::Product> {
  typedef AssemblableForward<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp,
                             L2Assemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>,
                             internal::L2AssemblableTraits, LocalEvaluation::Product> BaseType;
  template <class... Args>
  explicit L2Assemblable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};

template< class GridViewImp, class FieldImp, template<class, class, class> class OperatorTemplate>
class ProductForward
  : public ProductInterface< internal::DerivedType<GridViewImp, FieldImp, ProductForward<GridViewImp, FieldImp, OperatorTemplate > > >
{
public:
  typedef internal::DerivedType<GridViewImp, FieldImp,
                                ProductForward<GridViewImp, FieldImp, OperatorTemplate > > Traits;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const unsigned int                                   dimDomain = GridViewType::dimension;

  ProductForward(const GridViewType& grd_vw, const size_t over_integrate = 0)
    : grid_view_(grd_vw)
    , over_integrate_(over_integrate)
  {}

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template< class RR, int rRR, int rCR, class RS, int rRS, int rCS >
  FieldType apply2(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RR, rRR, rCR >& /*range*/,
                   const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RS, rRS, rCS >& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse< RR >::value), "Not implemented for this combination!");
  }

  template< int dimRangeRows, int dimRangeCols >
  FieldType apply2(const Stuff::LocalizableFunctionInterface
                      < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols >& range,
                   const Stuff::LocalizableFunctionInterface
                      < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols >& source) const
  {
    typedef Stuff::LocalizableFunctionInterface
        < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols > RangeType;
    OperatorTemplate< GridViewType, RangeType, RangeType >
        product_operator(grid_view_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class L2

template< class GridViewImp, class FieldImp = double>
using L2 = ProductForward<GridViewImp, FieldImp, L2Localizable>;

} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_HH
