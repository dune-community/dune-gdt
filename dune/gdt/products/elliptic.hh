// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ELLIPTIC_HH
#define DUNE_GDT_PRODUCTS_ELLIPTIC_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include "../localoperator/codim0.hh"
#include "../localevaluation/elliptic.hh"

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Products {


// forwards
template< class DiffusionType, class GridViewImp, class FieldImp = double >
class EllipticBase;

template< class DiffusionType, class GridViewImp, class RangeImp, class SourceImp = RangeImp, class FieldImp = double >
class EllipticLocalizable;

template< class DiffusionType
        , class MatrixImp
        , class RangeSpaceImp
        , class GridViewImp = typename RangeSpaceImp::GridViewType
        , class SourceSpaceImp = RangeSpaceImp >
class EllipticAssemblable;

template< class DiffusionType, class GridViewImp, class FieldImp = double >
class Elliptic;


namespace internal {


template< class DiffusionType, class GridViewImp, class FieldImp >
class EllipticBaseTraits
{
public:
  typedef GridViewImp GridViewType;
  typedef FieldImp    FieldType;
  typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionType > > LocalOperatorType;
}; // class EllipticBaseTraits


template< class DiffusionType, class GridViewImp, class FieldImp >
class EllipticBase
{
  typedef EllipticBaseTraits< DiffusionType, GridViewImp, FieldImp > Traits;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  EllipticBase(const DiffusionType& diffusion, const size_t over_integrate = 0)
    : diffusion_(diffusion)
    , local_operator_(over_integrate, diffusion_)
  {}

private:
  const DiffusionType& diffusion_;
protected:
  const LocalOperatorType local_operator_;
}; // class EllipticBase


template< class DiffusionType, class GridViewImp, class RangeImp, class SourceImp, class FieldImp >
class EllipticLocalizableTraits
  : public EllipticBaseTraits< DiffusionType, GridViewImp, FieldImp >
{
public:
  typedef EllipticLocalizable< DiffusionType, GridViewImp, RangeImp, SourceImp, FieldImp > derived_type;
  typedef RangeImp  RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp  FieldType;
};


template< class DiffusionType, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
class EllipticAssemblableTraits
  : public EllipticBaseTraits< DiffusionType, GridViewImp, typename MatrixImp::ScalarType >
{
public:
  typedef EllipticAssemblable< DiffusionType, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > derived_type;
  typedef RangeSpaceImp   RangeSpaceType;
  typedef SourceSpaceImp  SourceSpaceType;
  typedef MatrixImp       MatrixType;
}; // class EllipticAssemblableTraits


template< class DiffusionType, class GridViewImp, class FieldImp >
class EllipticTraits
{
public:
  typedef Elliptic< DiffusionType, GridViewImp, FieldImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
};


} // namespace internal


template< class DiffusionType, class GridViewImp, class RangeImp, class SourceImp, class FieldImp >
class EllipticLocalizable
  : public Products::LocalizableBase< internal::EllipticLocalizableTraits< DiffusionType, GridViewImp, RangeImp, SourceImp, FieldImp > >
  , public internal::EllipticBase< DiffusionType, GridViewImp, FieldImp >
{
  typedef Products::LocalizableBase
      < internal::EllipticLocalizableTraits< DiffusionType, GridViewImp, RangeImp, SourceImp, FieldImp > >
    LocalizableBaseType;
  typedef internal::EllipticBase< DiffusionType, GridViewImp, FieldImp > EllipticBaseType;
public:
  typedef internal::EllipticLocalizableTraits< DiffusionType, GridViewImp, RangeImp, SourceImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  EllipticLocalizable(const DiffusionType& diffusion,
                      const GridViewType& grd_vw,
                      const RangeType& rng,
                      const SourceType& src,
                      const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng, src)
    , EllipticBaseType(diffusion, over_integrate)
  {}

  EllipticLocalizable(const DiffusionType& diffusion,
                      const GridViewType& grd_vw,
                      const RangeType& rng,
                      const size_t over_integrate = 0)
    : LocalizableBaseType(grd_vw, rng)
    , EllipticBaseType(diffusion, over_integrate)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticLocalizable


template< class DiffusionType, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp >
class EllipticAssemblable
  : public Products::AssemblableBase< internal::EllipticAssemblableTraits< DiffusionType, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
  , public internal::EllipticBase< DiffusionType, GridViewImp, typename MatrixImp::ScalarType >
{
  typedef Products::AssemblableBase
      < internal::EllipticAssemblableTraits< DiffusionType, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp > >
    AssemblableBaseType;
  typedef internal::EllipticBase< DiffusionType, GridViewImp, typename MatrixImp::ScalarType > EllipticBaseType;
public:
  typedef internal::EllipticAssemblableTraits< DiffusionType, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp >
    Traits;
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

  EllipticAssemblable(const DiffusionType& diffusion,
                      MatrixType& matrix,
                      const RangeSpaceType& range_space,
                      const GridViewType& grid_view,
                      const SourceSpaceType& source_space)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , EllipticBaseType(diffusion)
  {}

  EllipticAssemblable(const DiffusionType& diffusion,
                      MatrixType& matrix,
                      const RangeSpaceType& range_space,
                      const GridViewType& grid_view)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , EllipticBaseType(diffusion)
  {}

  EllipticAssemblable(const DiffusionType& diffusion,
                      MatrixType& matrix,
                      const RangeSpaceType& range_space)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , EllipticBaseType(diffusion)
  {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticAssemblable


template< class DiffusionType, class GridViewImp, class FieldImp >
class Elliptic
  : public ProductInterface< internal::EllipticTraits< DiffusionType, GridViewImp, FieldImp > >
{
public:
  typedef internal::EllipticTraits< DiffusionType, GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  Elliptic(const DiffusionType& diffusion, const GridViewType& grd_vw)
    : diffusion_(diffusion)
    , grid_view_(grd_vw)
  {}

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template< class RR, int rRR, int rCR, class RS, int rRS, int rCS >
  FieldType apply2(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RR, rRR, rCR >& /*range*/,
                   const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RS, rRS, rCS >& /*source*/,
                   const size_t /*over_integrate*/ = 0) const
  {
    static_assert((Dune::AlwaysFalse< RR >::value), "Not implemented for this combination!");
  }

  template< int dimRangeRows, int dimRangeCols >
  FieldType apply2(const Stuff::LocalizableFunctionInterface
                      < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols >& range,
                   const Stuff::LocalizableFunctionInterface
                      < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols >& source,
                   const size_t over_integrate = 0) const
  {
    typedef Stuff::LocalizableFunctionInterface
        < EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols > FunctionType;
    EllipticLocalizable< DiffusionType, GridViewType, FunctionType, FunctionType, FieldType >
        product_operator(diffusion_, grid_view_, range, source, over_integrate);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const DiffusionType& diffusion_;
  const GridViewType& grid_view_;
}; // class Elliptic


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ELLIPTIC_HH
