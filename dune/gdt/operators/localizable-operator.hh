// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
#define DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/pattern.hh>
#include <dune/xt/functions/interfaces.hh>

#include <dune/gdt/local/assembler.hh>
#include <dune/gdt/assembler/wrapper.hh>
#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/type_traits.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GridLayerImp, class SourceImp, class RangeImp>
class LocalizableOperatorBase : public XT::Grid::Walker<GridLayerImp>
{
  typedef LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp> ThisType;
  typedef XT::Grid::Walker<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;

private:
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value, "Have to match!");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value, "Have to match!");
  static_assert(std::is_same<typename SourceType::DomainFieldType, typename GridLayerType::ctype>::value,
                "Have to match!");
  static_assert(std::is_same<typename RangeType::DomainFieldType, typename GridLayerType::ctype>::value,
                "Have to match!");
  static_assert(SourceType::dimDomain == GridLayerType::dimension, "Have to match!");
  static_assert(RangeType::dimDomain == GridLayerType::dimension, "Have to match!");

public:
  LocalizableOperatorBase(GridLayerType grd_layr, const SourceType& src, RangeType& rng)
    : BaseType(grd_layr)
    , source_(src)
    , range_(rng)
    , walked_(false)
  {
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

  using BaseType::grid_layer;
  using BaseType::append;

  template <class L>
  ThisType& append(
      const LocalOperatorInterface<L>& local_operator,
      const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where = new XT::Grid::ApplyOn::AllEntities<GridLayerType>())
  {
    typedef LocalOperatorApplicator<GridLayerType,
                                    typename LocalOperatorInterface<L>::derived_type,
                                    SourceType,
                                    RangeType>
        Applicator;
    local_operators_codim_0.emplace_back(
        new Applicator(grid_layer(), local_operator.as_imp(), source_, range_, *where));
    BaseType::append(*local_operators_codim_0.back(), where);
    return *this;
  } // ... append(...)

  template <class T>
  ThisType& append(const LocalCouplingOperatorInterface<T>& local_operator,
                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where =
                       new XT::Grid::ApplyOn::InnerIntersections<GridLayerType>(),
                   const XT::Common::Parameter& mu = {})
  {
    typedef LocalCouplingOperatorApplicator<GridLayerType,
                                            typename LocalCouplingOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType>
        Applicator;
    local_operators_codim_1.emplace_back(
        new Applicator(grid_layer(), local_operator.as_imp(), source_, range_, *where, mu));
    BaseType::append(*local_operators_codim_1.back(), where);
    return *this;
  } // ... append(...)

  template <class T>
  ThisType& append(const LocalCouplingOperatorInterface<T>& local_operator, const XT::Common::Parameter& mu)
  {
    return this->append(local_operator, new XT::Grid::ApplyOn::InnerIntersections<GridLayerType>(), mu);
  }

  template <class T>
  ThisType& append(const LocalBoundaryOperatorInterface<T>& local_operator,
                   const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where =
                       new XT::Grid::ApplyOn::BoundaryIntersections<GridLayerType>())
  {
    typedef LocalBoundaryOperatorApplicator<GridLayerType,
                                            typename LocalBoundaryOperatorInterface<T>::derived_type,
                                            SourceType,
                                            RangeType>
        Applicator;
    local_operators_codim_1.emplace_back(
        new Applicator(grid_layer(), local_operator.as_imp(), source_, range_, *where));
    BaseType::append(*local_operators_codim_1.back(), where);
    return *this;
  } // ... append(...)

  void apply(const bool use_tbb = false)
  {
    if (!walked_) {
      this->walk(use_tbb);
      walked_ = true;
    }
  } // ... apply(...)

protected:
  const SourceType& source_;
  RangeType& range_;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim0Object<GridLayerType>>> local_operators_codim_0;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim1Object<GridLayerType>>> local_operators_codim_1;
  bool walked_;
}; // class LocalizableOperatorBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
