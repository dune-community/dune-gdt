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

template< class MatrixImp
        , class RangeSpaceImp
        , class GridViewImp = typename RangeSpaceImp::GridViewType
        , class SourceSpaceImp = RangeSpaceImp >
class H1SemiAssemblable;


template< class GridViewImp, class FieldImp = double >
using H1SemiTraits = internal::L2BaseTraits<GridViewImp, FieldImp>;

template< class GridViewImp, class FieldImp >
using H1SemiBase = internal::L2Base<GridViewImp, FieldImp>;

template< class GridViewImp, class RangeImp, class SourceImp>
struct H1SemiLocalizable;

template< class GridViewImp, class RangeImp, class SourceImp, class DerivedImp >
using H1SemiLocalizableTraits = internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp, DerivedImp>;

template< class GridViewImp, class RangeImp, class SourceImp >
struct H1SemiLocalizable
    : public LocalizableForward<GridViewImp, RangeImp, SourceImp, H1SemiLocalizable<GridViewImp, RangeImp, SourceImp>,
                           H1SemiLocalizableTraits> {
  typedef LocalizableForward<GridViewImp, RangeImp, SourceImp, H1SemiLocalizable<GridViewImp, RangeImp, SourceImp>,
                        H1SemiLocalizableTraits> BaseType;
  template <class... Args>
  explicit H1SemiLocalizable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};

template< class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
class H1SemiAssemblableTraits
  : public H1SemiTraits< GridViewImp, typename MatrixImp::ScalarType >
{
public:
  typedef H1SemiAssemblable< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > derived_type;
  typedef RangeSpaceImp   RangeSpaceType;
  typedef SourceSpaceImp  SourceSpaceType;
  typedef MatrixImp       MatrixType;
private:
  friend class H1SemiAssemblable< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp >;
};


template< class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
class H1SemiAssemblable
  : public Products::AssemblableBase< H1SemiAssemblableTraits< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
  , public H1SemiBase< GridViewImp, typename MatrixImp::ScalarType >
{
  typedef Products::AssemblableBase< H1SemiAssemblableTraits< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
      AssemblableBaseType;
  typedef H1SemiBase< GridViewImp, typename MatrixImp::ScalarType > H1SemiBaseType;
public:
  typedef H1SemiAssemblableTraits< MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > Traits;
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

  H1SemiAssemblable(MatrixType& mtrx,
                    const RangeSpaceType& rng_spc,
                    const GridViewType& grd_vw,
                    const SourceSpaceType& src_spc,
                    const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_spc, grd_vw, src_spc)
    , H1SemiBaseType(over_integrate)
  {}

  H1SemiAssemblable(MatrixType& mtrx,
                    const RangeSpaceType& rng_spc,
                    const GridViewType& grd_vw,
                    const size_t over_integrate = 0)
    : AssemblableBaseType(mtrx, rng_spc, grd_vw, rng_spc)
    , H1SemiBaseType(over_integrate)
  {}

  H1SemiAssemblable(MatrixType& matrix, const RangeSpaceType& range_space, const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , H1SemiBaseType(over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class H1SemiAssemblable


template< class GridViewImp, class FieldImp = double >
class H1SemiGeneric;


template< class GridViewImp, class FieldImp = double >
class H1SemiGenericTraits
{
public:
  typedef H1SemiGeneric< GridViewImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
};


template< class GridViewImp, class FieldImp >
class H1SemiGeneric
  : public ProductInterface< H1SemiGenericTraits< GridViewImp, FieldImp > >
{
public:
  typedef H1SemiGenericTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  H1SemiGeneric(const GridViewType& grd_vw, const size_t over_integrate = 0)
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
        < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols > FunctionType;
    H1SemiLocalizable< GridViewType, FunctionType, FunctionType >
        product_operator(grid_view_, range, source, over_integrate_);
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
