// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PRODUCTS_HH
#define DUNE_GDT_OPERATOR_PRODUCTS_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/common/crtp.hh>

#include "../localoperator/codim0.hh"
#include "../localfunctional/codim0.hh"
#include "../localevaluation/product.hh"
#include "../localevaluation/elliptic.hh"

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace ProductOperator {


// forward, to be used in the traits
template< class GridPartImp, class FieldImp >
class L2Base;

template< class GridPartImp, class RangeImp, class SourceImp, class FieldImp = double >
class L2Localizable;

template< class GridPartImp, class RangeSpaceImp, class SourceSpaceImp, class MatrixImp >
class L2Assemblable;


template< class GridPartImp, class FieldImp = double >
class L2Traits
{
public:
  typedef GridPartImp GridPartType;
protected:
  typedef Stuff::Function::Constant<  typename GridPartType::template Codim< 0 >::EntityType,
                                      typename GridPartType::ctype,
                                      GridPartType::dimension, FieldImp, 1 >  FunctionType;
public:
  typedef LocalOperator::Codim0Integral< LocalEvaluation::Product< FunctionType > > LocalOperatorType;
private:
  friend class L2Base< GridPartImp, FieldImp >;
};

template< class GridPartImp, class FieldImp >
class L2Base
{
  typedef typename L2Traits< GridPartImp, FieldImp >::FunctionType       FunctionType;
  typedef typename L2Traits< GridPartImp, FieldImp >::LocalOperatorType  LocalOperatorType;

public:
  L2Base()
    : function_(1)
    , local_operator_(function_)
  {}

private:
  const FunctionType function_;
protected:
  const LocalOperatorType local_operator_;
}; // class L2Base


template< class GridPartImp, class RangeImp, class SourceImp, class FieldImp = double >
class L2LocalizableTraits
  : public L2Traits< GridPartImp, FieldImp >
{
public:
  typedef L2Localizable< GridPartImp, RangeImp, SourceImp, FieldImp > derived_type;
  typedef RangeImp    RangeType;
  typedef SourceImp   SourceType;
  typedef FieldImp    FieldType;
private:
  friend class L2Localizable< GridPartImp, RangeImp, SourceImp, FieldImp >;
};


template< class GridPartImp, class RangeImp, class SourceImp, class FieldImp >
class L2Localizable
  : public LocalizableProductBase< L2LocalizableTraits< GridPartImp, RangeImp, SourceImp, FieldImp > >
  , public L2Base< GridPartImp, FieldImp >
{
  typedef LocalizableProductBase< L2LocalizableTraits< GridPartImp, RangeImp, SourceImp, FieldImp > > BaseType;
public:
  typedef L2LocalizableTraits< GridPartImp, RangeImp, SourceImp, FieldImp > Traits;
  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  L2Localizable(const GridPartType& grid_part, const RangeType& range, const SourceType& source)
    : BaseType(grid_part, range, source)
  {}

  virtual const LocalOperatorType& local_operator() const DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Localizable


template< class GridPartImp, class RangeSpaceImp, class SourceSpaceImp, class MatrixImp >
class L2AssemblableTraits
  : public L2Traits< GridPartImp, typename MatrixImp::ScalarType >
{
public:
  typedef L2Assemblable< GridPartImp, RangeSpaceImp, SourceSpaceImp, MatrixImp > derived_type;
  typedef RangeSpaceImp   RangeSpaceType;
  typedef SourceSpaceImp  SourceSpaceType;
  typedef MatrixImp       MatrixType;
private:
  friend class L2Assemblable< GridPartImp, RangeSpaceImp, SourceSpaceImp, MatrixImp >;
};


template< class GridPartImp, class RangeSpaceImp, class SourceSpaceImp, class MatrixImp >
class L2Assemblable
  : public AssemblableProductBase< L2AssemblableTraits< GridPartImp, RangeSpaceImp, SourceSpaceImp, MatrixImp > >
  , public L2Base< GridPartImp, typename MatrixImp::ScalarType >
{
  typedef AssemblableProductBase< L2AssemblableTraits< GridPartImp, RangeSpaceImp, SourceSpaceImp, MatrixImp > >
    BaseType;
public:
  typedef L2AssemblableTraits< GridPartImp, RangeSpaceImp, SourceSpaceImp, MatrixImp > Traits;
  typedef typename Traits::GridPartType     GridPartType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::MatrixType       MatrixType;
private:
  typedef typename Traits::LocalOperatorType LocalOperatorType;

public:
  L2Assemblable(const GridPartType& grid_part, const RangeSpaceType& range_space, const SourceSpaceType& source_space)
    : BaseType(grid_part, range_space, source_space)
  {}

  virtual const LocalOperatorType& local_operator() const DS_FINAL
  {
    return this->local_operator_;
  }
}; // class L2Assemblable


template< class GridPartImp, class FieldImp = double >
class L2Generic;


template< class GridPartImp, class FieldImp = double >
class L2GenericTraits
{
public:
  typedef L2Generic< GridPartImp > derived_type;
  typedef GridPartImp GridPartType;
  typedef FieldImp FieldType;
};


template< class GridPartImp, class FieldImp >
class L2Generic
  : public ProductInterface< L2GenericTraits< GridPartImp, FieldImp > >
{
public:
  typedef L2GenericTraits< GridPartImp, FieldImp > Traits;
  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  L2Generic(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {}

  const GridPartType& grid_part() const
  {
    return grid_part_;
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
    L2Localizable< GridPartType, FunctionType, FunctionType, FieldType > product_operator(grid_part_, range, source);
    return product_operator.apply2();
  } // ... apply2(...)

private:
  const GridPartType& grid_part_;
}; // class L2Generic


template< class GridPartType >
class H1Semi
{
public:
  typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

  H1Semi(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {}

  template< class RR, int rRR, int rCR, class RS, int rRS, int rCS >
  void apply2(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RR, rRR, rCR >& /*range*/,
              const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, RS, rRS, rCS >& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse< RR >::value), "Not implemented for this combination!");
  }

  template< class RangeFieldType, int dimRangeRows, int dimRangeCols >
  RangeFieldType apply2(const Stuff::LocalizableFunctionInterface
                            < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols >&
                          range,
                        const Stuff::LocalizableFunctionInterface
                            < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRangeRows, dimRangeCols >&
                          source) const
  {
    // some preparations
    typedef Stuff::Function::Constant< EntityType, DomainFieldType, dimDomain, RangeFieldType, 1 > ConstantFunctionType;
    const ConstantFunctionType one(1);
    const LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< ConstantFunctionType > > local_operator(one);
    RangeFieldType ret = 0;
    DynamicMatrix< RangeFieldType > local_operator_result(1, 1, 0);
    std::vector< DynamicMatrix< RangeFieldType > > tmp_matrices(local_operator.numTmpObjectsRequired(),
                                                                DynamicMatrix< RangeFieldType >(1, 1, 0));

    // walk the grid
    const auto entity_it_end = grid_part_.template end< 0 >();
    for (auto entity_it = grid_part_.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto source_local_function = source.local_function(entity);
      const auto range_local_function = range.local_function(entity);
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
