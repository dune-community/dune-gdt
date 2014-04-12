// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_HH
#define DUNE_GDT_PRODUCTS_L2_HH

#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include "../localoperator/codim0.hh"
#include "../localevaluation/product.hh"

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Products {


// forward
template <class GridViewImp, class FunctionImp, class RangeImp, class SourceImp = RangeImp>
class WeightedL2Localizable;

// forward
template <class MatrixImp, class FunctionImp, class RangeSpaceImp,
          class GridViewImp = typename RangeSpaceImp::GridViewType, class SourceSpaceImp = RangeSpaceImp>
class WeightedL2Assemblable;

// forward
template <class GridViewImp, class FunctionImp>
class WeightedL2;

// forward
template <class GridViewImp, class RangeImp, class SourceImp = RangeImp>
class L2Localizable;

// forward
template <class MatrixImp, class RangeSpaceImp, class GridViewImp = typename RangeSpaceImp::GridViewType,
          class SourceSpaceImp                                    = RangeSpaceImp>
class L2Assemblable;

// forward
template <class GridViewImp, class FieldImp = double>
class L2;


/**
 *  \brief  Internal stuff that is usually not of interest to the user.
 */
namespace internal {


// forward
template <class GridViewImp, class FunctionImp>
class WeightedL2Base;

// forward
template <class GridViewImp, class FieldImp = double>
class L2Base;


template <class GridViewImp, class FunctionImp>
class WeightedL2BaseTraits
{
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, FunctionImp>::value,
                "FunctionImp has to be tagged as LocalizableFunction!");
  static_assert(std::is_base_of<GridView<typename GridViewImp::Traits>, GridViewImp>::value,
                "GridViewImp has to be derived from GridView!");
  static_assert(std::is_same<typename FunctionImp::DomainFieldType, typename GridViewImp::ctype>::value,
                "Types do not match!");
  static_assert(FunctionImp::dimDomain == GridViewImp::dimension, "Dimensions do not match!");

public:
  typedef GridViewImp GridViewType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Product<FunctionImp>> LocalOperatorType;
}; // class WeightedL2BaseTraits


template <class GridViewImp, class FunctionImp>
class WeightedL2Base
{
  typedef typename WeightedL2BaseTraits<GridViewImp, FunctionImp>::LocalOperatorType LocalOperatorType;

public:
  WeightedL2Base(const FunctionImp& function, const size_t over_integrate = 0)
    : local_operator_(over_integrate, function)
  {
  }

protected:
  const LocalOperatorType local_operator_;
}; // class WeightedL2Base


template <class GridViewImp, class FieldImp>
class L2BaseTraits
    : public WeightedL2BaseTraits<GridViewImp, Stuff::Function::Constant<
                                                   typename GridViewImp::template Codim<0>::Entity,
                                                   typename GridViewImp::ctype, GridViewImp::dimension, FieldImp, 1>>
{
  typedef Stuff::Function::Constant<typename GridViewImp::template Codim<0>::Entity, typename GridViewImp::ctype,
                                    GridViewImp::dimension, FieldImp, 1> FunctionType;

private:
  friend class L2Base<GridViewImp, FieldImp>;
}; // class L2BaseTraits


template <class GridViewImp, class FieldImp>
class L2Base
{
  typedef L2BaseTraits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::FunctionType FunctionType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  L2Base(const size_t over_integrate = 0)
    : function_(1)
    , local_operator_(over_integrate, function_)
  {
  }

protected:
  const FunctionType function_;
  const LocalOperatorType local_operator_;
}; // class L2Base


template <class GridViewImp, class FunctionImp, class RangeImp, class SourceImp>
class WeightedL2LocalizableTraits : public internal::WeightedL2BaseTraits<GridViewImp, FunctionImp>
{
public:
  typedef WeightedL2Localizable<GridViewImp, FunctionImp, RangeImp, SourceImp> derived_type;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef typename FunctionImp::RangeFieldType FieldType;
}; // class WeightedL2LocalizableTraits


template <class MatrixImp, class FunctionImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class WeightedL2AssemblableTraits : public internal::WeightedL2BaseTraits<GridViewImp, FunctionImp>
{
public:
  typedef WeightedL2Assemblable<MatrixImp, FunctionImp, RangeSpaceImp, GridViewImp, SourceSpaceImp> derived_type;
  typedef MatrixImp MatrixType;
  typedef FunctionImp FunctionType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
}; // class WeightedL2AssemblableTraits


template <class GridViewImp, class FunctionImp>
class WeightedL2Traits
{
public:
  typedef WeightedL2<GridViewImp, FunctionImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef typename FunctionImp::RangeFieldType FieldType;
}; // class WeightedL2Traits


template <class GridViewImp, class RangeImp, class SourceImp>
class L2LocalizableTraits : public internal::L2BaseTraits<GridViewImp, typename RangeImp::RangeFieldType>
{
  static_assert(std::is_same<typename RangeImp::RangeFieldType, typename SourceImp::RangeFieldType>::value,
                "Types do not match!");

public:
  typedef L2Localizable<GridViewImp, RangeImp, SourceImp> derived_type;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef typename RangeImp::RangeFieldType FieldType;
}; // class L2LocalizableTraits


template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class L2AssemblableTraits : public internal::L2BaseTraits<GridViewImp, typename RangeSpaceImp::RangeFieldType>
{
public:
  typedef L2Assemblable<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp> derived_type;
  typedef MatrixImp MatrixType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
}; // class L2AssemblableTraits


template <class GridViewImp, class FieldImp>
class L2Traits
{
public:
  typedef L2<GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class L2Traits


} // namespace internal


template <class GridViewImp, class FunctionImp, class RangeImp, class SourceImp>
class WeightedL2Localizable
    : public Products::LocalizableBase<internal::WeightedL2LocalizableTraits<GridViewImp, FunctionImp, RangeImp,
                                                                             SourceImp>>,
      public internal::WeightedL2Base<GridViewImp, FunctionImp>
{
  typedef Products::LocalizableBase<internal::WeightedL2LocalizableTraits<GridViewImp, FunctionImp, RangeImp,
                                                                          SourceImp>> LocalizableBaseType;
  typedef internal::WeightedL2Base<GridViewImp, FunctionImp> WeightedL2BaseType;

public:
  typedef internal::WeightedL2LocalizableTraits<GridViewImp, FunctionImp, RangeImp, SourceImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  WeightedL2Localizable(const GridViewType& grid_view, const FunctionImp& function, const RangeType& range,
                        const SourceType& source, const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range, source)
    , WeightedL2BaseType(function, over_integrate)
  {
  }

  WeightedL2Localizable(const GridViewType& grid_view, const FunctionImp& function, const RangeType& range,
                        const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range, range)
    , WeightedL2BaseType(function, over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class WeightedL2Localizable


template <class MatrixImp, class FunctionImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class WeightedL2Assemblable
    : public Products::AssemblableBase<internal::WeightedL2AssemblableTraits<MatrixImp, FunctionImp, RangeSpaceImp,
                                                                             GridViewImp, SourceSpaceImp>>,
      public internal::WeightedL2Base<GridViewImp, FunctionImp>
{
  typedef Products::AssemblableBase<internal::WeightedL2AssemblableTraits<MatrixImp, FunctionImp, RangeSpaceImp,
                                                                          GridViewImp, SourceSpaceImp>>
      AssemblableBaseType;
  typedef internal::WeightedL2Base<GridViewImp, FunctionImp> WeightedL2BaseType;

public:
  typedef internal::WeightedL2AssemblableTraits<MatrixImp, FunctionImp, RangeSpaceImp, GridViewImp, SourceSpaceImp>
      Traits;
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

  WeightedL2Assemblable(MatrixType& matrix, const FunctionImp& function, const RangeSpaceType& range_space,
                        const GridViewType& grid_view, const SourceSpaceType& source_space,
                        const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , WeightedL2BaseType(function, over_integrate)
  {
  }

  WeightedL2Assemblable(MatrixType& matrix, const FunctionImp& function, const RangeSpaceType& range_space,
                        const GridViewType& grid_view, const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , WeightedL2BaseType(function, over_integrate)
  {
  }

  WeightedL2Assemblable(MatrixType& matrix, const FunctionImp& function, const RangeSpaceType& range_space,
                        const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , WeightedL2BaseType(function, over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class WeightedL2Assemblable


template <class GridViewImp, class FunctionImp>
class WeightedL2 : public ProductInterface<internal::WeightedL2Traits<GridViewImp, FunctionImp>>
{
public:
  typedef internal::WeightedL2Traits<GridViewImp, FunctionImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  WeightedL2(const GridViewType& grid_view, const FunctionImp& function, const size_t over_integrate = 0)
    : grid_view_(grid_view)
    , function_(function)
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
            RangeType;
    WeightedL2Localizable<GridViewType, FunctionImp, RangeType, RangeType> product_operator(
        grid_view_, function_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const FunctionImp& function_;
  const size_t over_integrate_;
}; // class WeightedL2


template <class GridViewImp, class RangeImp, class SourceImp>
class L2Localizable : public Products::LocalizableBase<internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp>>,
                      public internal::L2Base<GridViewImp, typename RangeImp::RangeFieldType>
{
  typedef Products::LocalizableBase<internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp>>
      LocalizableBaseType;
  typedef internal::L2Base<GridViewImp, typename RangeImp::RangeFieldType> L2BaseType;

public:
  typedef internal::L2LocalizableTraits<GridViewImp, RangeImp, SourceImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;

private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  L2Localizable(const GridViewType& grid_view, const RangeType& range, const SourceType& source,
                const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range, source)
    , L2BaseType(over_integrate)
  {
  }

  L2Localizable(const GridViewType& grid_view, const RangeType& range, const size_t over_integrate = 0)
    : LocalizableBaseType(grid_view, range, range)
    , L2BaseType(over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Localizable


template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp>
class L2Assemblable : public Products::AssemblableBase<internal::L2AssemblableTraits<MatrixImp, RangeSpaceImp,
                                                                                     GridViewImp, SourceSpaceImp>>,
                      public internal::L2Base<GridViewImp, typename RangeSpaceImp::RangeFieldType>
{
  typedef Products::AssemblableBase<internal::L2AssemblableTraits<MatrixImp, RangeSpaceImp, GridViewImp,
                                                                  SourceSpaceImp>> AssemblableBaseType;
  typedef internal::L2Base<GridViewImp, typename RangeSpaceImp::RangeFieldType> L2BaseType;

public:
  typedef internal::L2AssemblableTraits<MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp> Traits;
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

  L2Assemblable(MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view,
                const SourceSpaceType& source_space, const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, source_space)
    , L2BaseType(over_integrate)
  {
  }

  L2Assemblable(MatrixType& matrix, const RangeSpaceType& range_space, const GridViewType& grid_view,
                const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, grid_view, range_space)
    , L2BaseType(over_integrate)
  {
  }

  L2Assemblable(MatrixType& matrix, const RangeSpaceType& range_space, const size_t over_integrate = 0)
    : AssemblableBaseType(matrix, range_space, *(range_space.grid_view()), range_space)
    , L2BaseType(over_integrate)
  {
  }

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Assemblable


template <class GridViewImp, class FieldImp>
class L2 : public ProductInterface<internal::L2Traits<GridViewImp, FieldImp>>
{
public:
  typedef internal::L2Traits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  L2(const GridViewType& grid_view, const size_t over_integrate = 0)
    : grid_view_(grid_view)
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
            RangeType;
    L2Localizable<GridViewType, RangeType, RangeType> product_operator(grid_view_, range, source, over_integrate_);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class L2


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_HH
