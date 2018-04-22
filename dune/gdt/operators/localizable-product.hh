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

#ifndef DUNE_GDT_OPERATORS_LOCALIZABLE_PRODUCT_HH
#define DUNE_GDT_OPERATORS_LOCALIZABLE_PRODUCT_HH

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


/**
 * \todo Check parallel case: there is probably/definitely communication missing in apply2!
 */
template <class GridLayerImp,
          class RangeImp,
          class SourceImp = RangeImp,
          class FieldImp = typename RangeImp::RangeFieldType>
class LocalizableProductBase : public XT::Grid::Walker<GridLayerImp>
{
  typedef LocalizableProductBase<GridLayerImp, RangeImp, SourceImp> ThisType;
  typedef XT::Grid::Walker<GridLayerImp> BaseType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;
  typedef typename RangeType::LocalfunctionType LocalRangeType;
  typedef typename SourceType::LocalfunctionType LocalSourceType;

private:
  static_assert(XT::Functions::is_localizable_function<SourceType>::value, "");
  static_assert(XT::Functions::is_localizable_function<RangeType>::value, "");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value, "");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value, "");
  static_assert(std::is_same<typename SourceType::DomainFieldType, typename GridLayerType::ctype>::value, "");
  static_assert(std::is_same<typename RangeType::DomainFieldType, typename GridLayerType::ctype>::value, "");
  static_assert(SourceType::dimDomain == GridLayerType::dimension, "");
  static_assert(RangeType::dimDomain == GridLayerType::dimension, "");

public:
  LocalizableProductBase(GridLayerType grd_layr,
                         const RangeType& rng,
                         const SourceType& src,
                         const XT::Common::Parameter& param = {})
    : BaseType(grd_layr)
    , range_(rng)
    , source_(src)
    , result_(0.)
    , walked_(false)
    , param_(param)
  {
  }

  LocalizableProductBase(GridLayerType grd_layr, const RangeType& rng, const XT::Common::Parameter& param = {})
    : BaseType(grd_layr)
    , range_(rng)
    , source_(rng)
    , result_(0.)
    , walked_(false)
    , param_(param)
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

  using BaseType::grid_layer;
  using BaseType::append;

  ThisType& append(const LocalVolumeTwoFormInterface<LocalRangeType, LocalSourceType, FieldType>& local_volume_twoform,
                   const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where =
                       new XT::Grid::ApplyOn::PartitionSetEntities<GridLayerType, Partitions::InteriorBorder>())
  {
    local_volume_twoforms_.emplace_back(
        new LocalVolumeTwoFormAccumulatorFunctor<GridLayerType, RangeType, SourceType, FieldType>(
            this->grid_layer_, local_volume_twoform, range_, source_, result_, where->copy()));
    BaseType::append(*local_volume_twoforms_.back(), where);
    return *this;
  } // ... append(...)

  FieldType compute_locally(const EntityType& entity) const
  {
    FieldType local_result = 0.;
    for (const auto& local_volume_twoform : local_volume_twoforms_)
      local_result += local_volume_twoform->compute_locally(entity);
    return local_result;
  }

  FieldType apply2()
  {
    if (!walked_) {
      this->walk();
      walked_ = true;
    }
    return result_;
  }

  FieldType result()
  {
    return result_;
  }

  const XT::Common::Parameter& parameter() const
  {
    return param_;
  }

protected:
  const RangeType& range_;
  const SourceType& source_;
  FieldType result_;
  std::vector<std::unique_ptr<XT::Grid::internal::Codim0ReturnObject<GridLayerType, FieldType>>> local_volume_twoforms_;
  bool walked_;
  const XT::Common::Parameter param_;
}; // class LocalizableProductBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LOCALIZABLE_PRODUCT_HH
