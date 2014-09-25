// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_WEIGHTED_HH
#define DUNE_GDT_PRODUCTS_L2_WEIGHTED_HH

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

template< class GridViewImp, class FunctionImp, class RangeImp, class SourceImp >
class WeightedL2Localizable
  : public Products::LocalizableBase< internal::WeightedL2LocalizableTraits< GridViewImp, FunctionImp, RangeImp, SourceImp > >
  , public internal::WeightedL2Base< GridViewImp, FunctionImp >
{
  typedef Products::LocalizableBase< internal::WeightedL2LocalizableTraits< GridViewImp, FunctionImp, RangeImp, SourceImp > >
    LocalizableBaseType;
  typedef internal::WeightedL2Base< GridViewImp, FunctionImp > WeightedL2BaseType;
public:
  typedef internal::WeightedL2LocalizableTraits< GridViewImp, FunctionImp, RangeImp, SourceImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  WeightedL2Localizable(const GridViewType& grd_vw,
                        const FunctionImp& function,
                        const RangeType& rng,
                        const SourceType& src,
                        const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, src)
    , WeightedL2BaseType(function, over_integrate)
  {}

  WeightedL2Localizable(const GridViewType& grd_vw,
                        const FunctionImp& function,
                        const RangeType& rng,
                        const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, rng)
    , WeightedL2BaseType(function, over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class WeightedL2Localizable


template< class MatrixImp, class FunctionImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
class WeightedL2Assemblable
  : public Products::AssemblableBase< internal::WeightedL2AssemblableTraits< MatrixImp, FunctionImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
  , public internal::WeightedL2Base< GridViewImp, FunctionImp >
{
  typedef Products::AssemblableBase< internal::WeightedL2AssemblableTraits< MatrixImp, FunctionImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
    AssemblableBaseType;
  typedef internal::WeightedL2Base< GridViewImp, FunctionImp > WeightedL2BaseType;
public:
  typedef internal::WeightedL2AssemblableTraits< MatrixImp, FunctionImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > Traits;
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

  WeightedL2Assemblable(MatrixType& matrix,
                        const FunctionImp& function,
                        const RangeSpaceType& range_space,
                        const GridViewType& grid_view,
                        const SourceSpaceType& source_space,
                        const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , WeightedL2BaseType(function, over_integrate)
  {}

  WeightedL2Assemblable(MatrixType& matrix,
                        const FunctionImp& function,
                        const RangeSpaceType& range_space,
                        const GridViewType& grid_view,
                        const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , WeightedL2BaseType(function, over_integrate)
  {}

  WeightedL2Assemblable(MatrixType& matrix,
                        const FunctionImp& function,
                        const RangeSpaceType& range_space,
                        const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , WeightedL2BaseType(function, over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class WeightedL2Assemblable


template< class GridViewImp, class FunctionImp >
class WeightedL2
  : public ProductInterface< internal::WeightedL2Traits< GridViewImp, FunctionImp > >
{
public:
  typedef internal::WeightedL2Traits< GridViewImp, FunctionImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const unsigned int                                   dimDomain = GridViewType::dimension;

  WeightedL2(const GridViewType& grd_vw, const FunctionImp& function, const size_t over_integrate = 0)
    : grid_view_(grd_vw)
    , function_(function)
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
    WeightedL2Localizable< GridViewType, FunctionImp, RangeType, RangeType >
        product_operator(grid_view_, function_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const FunctionImp& function_;
  const size_t over_integrate_;
}; // class WeightedL2


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_WEIGHTED_HH
