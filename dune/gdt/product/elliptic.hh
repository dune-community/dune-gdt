// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCT_ELLIPTIC_HH
#define DUNE_GDT_PRODUCT_ELLIPTIC_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include "../localoperator/codim0.hh"
#include "../localevaluation/elliptic.hh"

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Product {


// forward, to be used in the traits
template <class DiffusionImp, class GridViewImp, class FieldImp = double>
class EllipticBase;


template <class DiffusionImp, class GridViewImp, class RangeImp, class SourceImp = RangeImp, class FieldImp = double>
class EllipticLocalizable;


template <class DiffusionImp, class MatrixImp, class RangeSpaceImp,
          class GridViewImp = typename RangeSpaceImp::GridViewType, class SourceSpaceImp = RangeSpaceImp>
class EllipticAssemblable;


template <class DiffusionImp, class GridViewImp, class FieldImp>
class EllipticBaseTraits
{
public:
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;

protected:
  typedef DiffusionImp DiffusionType;

public:
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> LocalOperatorType;

private:
  friend class EllipticBase<DiffusionImp, GridViewImp, FieldImp>;
};


template <class DiffusionImp, class GridViewImp, class FieldImp>
class EllipticBase
{
  typedef typename EllipticBaseTraits<DiffusionImp, GridViewImp, FieldImp>::DiffusionType DiffusionType;
  typedef typename EllipticBaseTraits<DiffusionImp, GridViewImp, FieldImp>::LocalOperatorType LocalOperatorType;

public:
  EllipticBase(const DiffusionType& diffusion, const size_t over_integrate = 0)
    : diffusion_(diffusion)
    , local_operator_(over_integrate, diffusion_)
  {
  }

private:
  const DiffusionType& diffusion_;

protected:
  const LocalOperatorType local_operator_;
}; // class EllipticBase


template <class DiffusionImp, class GridViewImp, class RangeImp, class SourceImp, class FieldImp>
class EllipticLocalizableTraits : public EllipticBaseTraits<DiffusionImp, GridViewImp, FieldImp>
{
public:
  typedef EllipticLocalizable<DiffusionImp, GridViewImp, RangeImp, SourceImp, FieldImp> derived_type;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;

private:
  friend class EllipticLocalizable<DiffusionImp, GridViewImp, RangeImp, SourceImp, FieldImp>;
};


template <class DiffusionImp, class GridViewImp, class RangeImp, class SourceImp, class FieldImp>
class EllipticLocalizable : public Product::LocalizableBase<EllipticLocalizableTraits<DiffusionImp, GridViewImp,
                                                                                      RangeImp, SourceImp, FieldImp>>,
                            public EllipticBase<DiffusionImp, GridViewImp, FieldImp>
{
  typedef Product::LocalizableBase<EllipticLocalizableTraits<DiffusionImp, GridViewImp, RangeImp, SourceImp, FieldImp>>
      LocalizableBaseType;
  typedef EllipticBase<DiffusionImp, GridViewImp, FieldImp> EllipticBaseType;

public:
  typedef EllipticLocalizableTraits<DiffusionImp, GridViewImp, RangeImp, SourceImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  EllipticLocalizable(const DiffusionImp& diffusion, const GridViewType& grid_view, const RangeType& range,
                      const SourceType& source, const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range, source)
    , EllipticBaseType(diffusion, over_integrate)
  {
  }

  EllipticLocalizable(const DiffusionImp& diffusion, const GridViewType& grid_view, const RangeType& range,
                      const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range)
    , EllipticBaseType(diffusion, over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticLocalizable


template <class DiffusionImp, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class EllipticAssemblableTraits : public EllipticBaseTraits<DiffusionImp, GridViewImp, typename MatrixImp::ScalarType>
{
public:
  typedef EllipticAssemblable<DiffusionImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp> derived_type;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef MatrixImp MatrixType;

private:
  friend class EllipticAssemblable<DiffusionImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>;
};


template <class DiffusionImp, class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class EllipticAssemblable
    : public Product::AssemblableBase<EllipticAssemblableTraits<DiffusionImp, MatrixImp, RangeSpaceImp, GridViewImp,
                                                                SourceSpaceImp>>,
      public EllipticBase<DiffusionImp, GridViewImp, typename MatrixImp::ScalarType>
{
  typedef Product::AssemblableBase<EllipticAssemblableTraits<DiffusionImp, MatrixImp, RangeSpaceImp, GridViewImp,
                                                             SourceSpaceImp>> AssemblableBaseType;
  typedef EllipticBase<DiffusionImp, GridViewImp, typename MatrixImp::ScalarType> EllipticBaseType;

public:
  typedef EllipticAssemblableTraits<DiffusionImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::MatrixType MatrixType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  using AssemblableBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  EllipticAssemblable(const DiffusionImp& diffusion, MatrixType& matrix, const RangeSpaceType& range_space,
                      const GridViewType& grid_view, const SourceSpaceType& source_space)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , EllipticBaseType(diffusion)
  {
  }

  EllipticAssemblable(const DiffusionImp& diffusion, MatrixType& matrix, const RangeSpaceType& range_space,
                      const GridViewType& grid_view)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , EllipticBaseType(diffusion)
  {
  }

  EllipticAssemblable(const DiffusionImp& diffusion, MatrixType& matrix, const RangeSpaceType& range_space)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , EllipticBaseType(diffusion)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class EllipticAssemblable


template <class DiffusionImp, class GridViewImp, class FieldImp = double>
class Elliptic;


template <class DiffusionImp, class GridViewImp, class FieldImp>
class EllipticTraits
{
public:
  typedef Elliptic<DiffusionImp, GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
};


template <class DiffusionImp, class GridViewImp, class FieldImp>
class Elliptic : public ProductInterface<EllipticTraits<DiffusionImp, GridViewImp, FieldImp>>
{
public:
  typedef EllipticTraits<DiffusionImp, GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  Elliptic(const DiffusionImp& diffusion, const GridViewType& grid_view)
    : diffusion_(diffusion)
    , grid_view_(grid_view)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/,
         const size_t /*over_integrate*/ = 0) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <int dimRangeRows, int dimRangeCols>
  FieldType apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& range,
                   const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType,
                                                             dimRangeRows, dimRangeCols>& source,
                   const size_t over_integrate = 0) const
  {
    typedef Stuff::
        LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, dimRangeRows, dimRangeCols>
            FunctionType;
    EllipticLocalizable<DiffusionImp, GridViewType, FunctionType, FunctionType, FieldType> product_operator(
        diffusion_, grid_view_, range, source, over_integrate);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const DiffusionImp& diffusion_;
  const GridViewType& grid_view_;
}; // class Elliptic


} // namespace Operator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCT_ELLIPTIC_HH
