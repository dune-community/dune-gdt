// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PRODUCTS_HH
#define DUNE_GDT_OPERATOR_PRODUCTS_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "../localoperator/codim0.hh"
#include "../localfunctional/codim0.hh"
#include "../localevaluation/product.hh"
#include "../localevaluation/elliptic.hh"

namespace Dune {
namespace GDT {
namespace ProductOperator {


template <class GridPartType>
class L2
{
public:
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  L2(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  void apply2(
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <class RangeFieldType, int dimRangeRows, int dimRangeCols>
  RangeFieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows,
                                                   dimRangeCols>& range,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows,
                                                   dimRangeCols>& source) const
  {
    typedef Stuff::
        LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols>
            RangeType;
    // some preparations
    const LocalFunctional::Codim0Integral<LocalEvaluation::Product<RangeType>> local_functional(range);
    RangeFieldType ret = 0;
    DynamicVector<RangeFieldType> local_functional_result(1, 0);
    std::vector<DynamicVector<RangeFieldType>> tmp_vectors(local_functional.numTmpObjectsRequired(),
                                                           DynamicVector<RangeFieldType>(1, 0));

    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local function
      const auto source_local_function = source.local_function(entity);
      // apply local functional
      local_functional.apply(*source_local_function, local_functional_result, tmp_vectors);
      assert(local_functional_result.size() == 1);
      ret += local_functional_result[0];
    } // walk the grid
    return ret;
  } // ... apply2(...)

private:
  const GridPartType& grid_part_;
}; // class L2


template <class GridPartType>
class H1Semi
{
public:
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  H1Semi(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  void apply2(
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  template <class RangeFieldType, int dimRangeRows, int dimRangeCols>
  RangeFieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows,
                                                   dimRangeCols>& range,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows,
                                                   dimRangeCols>& source) const
  {
    // some preparations
    typedef Stuff::Function::Constant<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1> ConstantFunctionType;
    const ConstantFunctionType one(1);
    const LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<ConstantFunctionType>> local_operator(one);
    RangeFieldType ret = 0;
    DynamicMatrix<RangeFieldType> local_operator_result(1, 1, 0);
    std::vector<DynamicMatrix<RangeFieldType>> tmp_matrices(local_operator.numTmpObjectsRequired(),
                                                            DynamicMatrix<RangeFieldType>(1, 1, 0));

    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto source_local_function = source.local_function(entity);
      const auto range_local_function  = range.local_function(entity);
      // apply local operator
      local_operator.apply(*source_local_function, *range_local_function, local_operator_result, tmp_matrices);
      assert(local_operator_result.rows() == 1);
      assert(local_operator_result.cols() == 1);
      ret += local_operator_result[0][0];
    } // walk the grid
    return ret;
  } // ... apply2(...)

private:
  const GridPartType& grid_part_;
}; // class H1Semi


} // namespace ProductOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PRODUCTS_HH
