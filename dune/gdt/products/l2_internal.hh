// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_L2_BASE_HH
#define DUNE_GDT_PRODUCTS_L2_BASE_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/grid/common/gridview.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/products/base.hh>

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

//! Internal stuff that is usually not of interest to the user.
namespace internal {

// forward
template <class GridViewImp, class FunctionImp>
class WeightedL2Base;

// forward
template <class GridViewImp, class FieldImp, template <class> class LocalEvaluationTemplate>
class L2Base;


template <class GridViewImp, class FunctionImp,
          template <class> class LocalEvaluationTemplate = LocalEvaluation::Product>
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
  typedef LocalOperator::Codim0Integral<LocalEvaluationTemplate<FunctionImp>> LocalOperatorType;
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


template <class GridViewImp, class FieldImp, template <class> class LocalEvaluationTemplate>
struct L2BaseTraits
    : public WeightedL2BaseTraits<GridViewImp, Stuff::Functions::Constant<
                                                   typename GridViewImp::template Codim<0>::Entity,
                                                   typename GridViewImp::ctype, GridViewImp::dimension, FieldImp, 1>,
                                  LocalEvaluationTemplate>
{
  typedef Stuff::Functions::Constant<typename GridViewImp::template Codim<0>::Entity, typename GridViewImp::ctype,
                                     GridViewImp::dimension, FieldImp, 1> FunctionType;
}; // class L2BaseTraits


template <class GridViewImp, class FieldImp, template <class> class LocalEvaluationTemplate>
class L2Base
{
  typedef L2BaseTraits<GridViewImp, FieldImp, LocalEvaluationTemplate> Traits;
  typedef typename Traits::FunctionType FunctionType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  L2Base(const size_t over_integrate = 0)
    : function_(1)
    , local_operator_(over_integrate, function_)
  {
  }

private:
  const FunctionType function_;

protected:
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


template <class GridViewImp, class RangeImp, class SourceImp, class DerivedImp,
          template <class> class LocalEvaluationTemplate>
class L2LocalizableTraits
    : public internal::L2BaseTraits<GridViewImp, typename RangeImp::RangeFieldType, LocalEvaluationTemplate>
{
  static_assert(std::is_same<typename RangeImp::RangeFieldType, typename SourceImp::RangeFieldType>::value,
                "Types do not match!");

public:
  typedef DerivedImp derived_type;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef typename RangeImp::RangeFieldType FieldType;
}; // class L2LocalizableTraits


template <class MatrixImp, class RangeSpaceImp, class GridViewImp, class SourceSpaceImp, class DerivedImp,
          template <class> class LocalEvaluationTemplate = LocalEvaluation::Product>
struct L2AssemblableTraits
    : public internal::L2BaseTraits<GridViewImp, typename RangeSpaceImp::RangeFieldType, LocalEvaluationTemplate>
{
  typedef DerivedImp derived_type;
  typedef MatrixImp MatrixType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
}; // class L2AssemblableTraits


template <class GridViewImp, class FieldImp, class DerivedImp>
struct DerivedType
{
public:
  typedef DerivedImp derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class L2Traits

} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_L2_BASE_HH
