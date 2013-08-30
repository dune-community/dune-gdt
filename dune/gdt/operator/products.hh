// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PRODUCTS_HH
#define DUNE_GDT_OPERATOR_PRODUCTS_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "../basefunctionset/wrapper.hh"
#include "../localfunctional/codim0.hh"
#include "../localevaluation/product.hh"
#include "../localoperator/codim0.hh"
#include "../localevaluation/elliptic.hh"

namespace Dune {
namespace GDT {
namespace ProductOperator {


template <class GridPartType>
class L2
{
public:
  L2(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class RangeType, class SourceType>
  typename RangeType::template LocalFunction<typename GridPartType::template Codim<0>::EntityType>::Type::RangeFieldType
  apply2(const RangeType& range, const SourceType& source) const
  {
    // checks
    static_assert(std::is_base_of<Stuff::LocalizableFunction, SourceType>::value,
                  "SourceType has to be derived from Stuff::LocalizableFunction!");
    static_assert(std::is_base_of<Stuff::LocalizableFunction, RangeType>::value,
                  "RangeType has to be derived from Stuff::LocalizableFunction!");
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
    typedef typename SourceType::template LocalFunction<EntityType>::Type SourceLocalFunctionType;
    typedef typename RangeType::template LocalFunction<EntityType>::Type RangeLocalFunctionType;
    typedef typename SourceLocalFunctionType::DomainFieldType DomainFieldType;
    static_assert(std::is_same<DomainFieldType, typename RangeLocalFunctionType::DomainFieldType>::value,
                  "DomainFieldType of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimDomain = SourceLocalFunctionType::dimDomain;
    static_assert(dimDomain == RangeLocalFunctionType::dimDomain,
                  "dimDomain of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    typedef typename SourceLocalFunctionType::RangeFieldType RangeFieldType;
    static_assert(std::is_same<RangeFieldType, typename RangeLocalFunctionType::RangeFieldType>::value,
                  "RangeFieldType of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimRangeRows = SourceLocalFunctionType::dimRangeRows;
    static_assert(dimRangeRows == RangeLocalFunctionType::dimRangeRows,
                  "dimRangeRows of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimRangeCols = SourceLocalFunctionType::dimRangeCols;
    static_assert(dimRangeCols == RangeLocalFunctionType::dimRangeCols,
                  "dimRangeCols of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    // some preparations
    typedef LocalFunctional::Codim0Integral<LocalEvaluation::Product<RangeType>> LocalFunctionalType;
    const LocalFunctionalType local_functional(range);
    RangeFieldType ret = 0;
    DynamicVector<RangeFieldType> local_functional_result(1, 0);
    std::vector<DynamicVector<RangeFieldType>> tmp_vectors(local_functional.numTmpObjectsRequired(),
                                                           DynamicVector<RangeFieldType>(1, 0));

    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local function
      const auto source_local_function = source.localFunction(entity);
      typedef BaseFunctionSet::LocalFunctionWrapper<SourceLocalFunctionType,
                                                    DomainFieldType,
                                                    dimDomain,
                                                    RangeFieldType,
                                                    dimRangeRows,
                                                    dimRangeCols> LocalSourceWrapperType;
      const LocalSourceWrapperType local_source_wrapper(source_local_function);
      // apply local functional
      local_functional.apply(local_source_wrapper, local_functional_result, tmp_vectors);
      assert(local_functional_result.size() == 1);
      ret += local_functional_result[0];
    } // walk the grid
    return ret;
  }

private:
  const GridPartType& grid_part_;
}; // class L2


template <class GridPartType>
class H1
{
public:
  H1(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class RangeType, class SourceType>
  typename RangeType::template LocalFunction<typename GridPartType::template Codim<0>::EntityType>::Type::RangeFieldType
  apply2(const RangeType& range, const SourceType& source) const
  {
    // checks
    static_assert(std::is_base_of<Stuff::LocalizableFunction, SourceType>::value,
                  "SourceType has to be derived from Stuff::LocalizableFunction!");
    static_assert(std::is_base_of<Stuff::LocalizableFunction, RangeType>::value,
                  "RangeType has to be derived from Stuff::LocalizableFunction!");
    typedef typename GridPartType::template Codim<0>::EntityType EntityType;
    typedef typename SourceType::template LocalFunction<EntityType>::Type SourceLocalFunctionType;
    typedef typename RangeType::template LocalFunction<EntityType>::Type RangeLocalFunctionType;
    typedef typename SourceLocalFunctionType::DomainFieldType DomainFieldType;
    static_assert(std::is_same<DomainFieldType, typename RangeLocalFunctionType::DomainFieldType>::value,
                  "DomainFieldType of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimDomain = SourceLocalFunctionType::dimDomain;
    static_assert(dimDomain == RangeLocalFunctionType::dimDomain,
                  "dimDomain of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    typedef typename SourceLocalFunctionType::RangeFieldType RangeFieldType;
    static_assert(std::is_same<RangeFieldType, typename RangeLocalFunctionType::RangeFieldType>::value,
                  "RangeFieldType of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimRangeRows = SourceLocalFunctionType::dimRangeRows;
    static_assert(dimRangeRows == RangeLocalFunctionType::dimRangeRows,
                  "dimRangeRows of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    static const unsigned int dimRangeCols = SourceLocalFunctionType::dimRangeCols;
    static_assert(dimRangeCols == RangeLocalFunctionType::dimRangeCols,
                  "dimRangeCols of SourceLocalFunctionType and RangeLocalFunctionType do not match!");
    // some preparations
    typedef Stuff::FunctionConstant<DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
        ConstantFunctionType;
    const ConstantFunctionType one(1);
    typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<ConstantFunctionType>> LocalOperatorType;
    const LocalOperatorType local_operator(one);
    RangeFieldType ret = 0;
    DynamicMatrix<RangeFieldType> local_operator_result(1, 1, 0);
    std::vector<DynamicMatrix<RangeFieldType>> tmp_matrices(local_operator.numTmpObjectsRequired(),
                                                            DynamicMatrix<RangeFieldType>(1, 1, 0));

    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto source_local_function = source.localFunction(entity);
      typedef BaseFunctionSet::LocalFunctionWrapper<SourceLocalFunctionType,
                                                    DomainFieldType,
                                                    dimDomain,
                                                    RangeFieldType,
                                                    dimRangeRows,
                                                    dimRangeCols> LocalSourceWrapperType;
      const LocalSourceWrapperType local_source_wrapper(source_local_function);
      const auto range_local_function = range.localFunction(entity);
      typedef BaseFunctionSet::LocalFunctionWrapper<RangeLocalFunctionType,
                                                    DomainFieldType,
                                                    dimDomain,
                                                    RangeFieldType,
                                                    dimRangeRows,
                                                    dimRangeCols> LocalRangeWrapperType;
      const LocalRangeWrapperType local_range_wrapper(range_local_function);
      // apply local operator
      local_operator.apply(local_source_wrapper, local_range_wrapper, local_operator_result, tmp_matrices);
      assert(local_operator_result.rows() == 1);
      assert(local_operator_result.cols() == 1);
      ret += local_operator_result[0][0];
    } // walk the grid
    return ret;
  }

private:
  const GridPartType& grid_part_;
}; // class H1


} // namespace ProductOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PRODUCTS_HH
