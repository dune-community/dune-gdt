// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROJECTIONS_L2_LOCAL_HH
#define DUNE_GDT_PROJECTIONS_L2_LOCAL_HH

#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/l2-projection.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/spaces/dg/interface.hh>

namespace Dune {
namespace GDT {


// forward
template <class GridLayerImp, class FieldImp = double>
class L2LocalProjectionOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class L2LocalProjectionOperatorTraits
{
public:
  typedef L2LocalProjectionOperator<GridLayerImp, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class SourceImp, class RangeImp>
class L2LocalProjectionLocalizableOperator : public LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp>
{
  typedef LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp> BaseType;
  typedef LocalL2ProjectionOperator LocalOperatorType;

public:
  using typename BaseType::GridLayerType;
  using typename BaseType::SourceType;
  using typename BaseType::RangeType;

  template <class... Args>
  explicit L2LocalProjectionLocalizableOperator(const size_t over_integrate, Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_operator_(over_integrate)
  {
    this->append(local_operator_);
    issue_warning(this->range().space());
  }

  template <class... Args>
  explicit L2LocalProjectionLocalizableOperator(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
    , local_operator_()
  {
    this->append(local_operator_);
    issue_warning(this->range().space());
  }

private:
  template <class S, bool warn>
  struct Warning
  {
    static inline void issue()
    {
    }
  };

  template <class S>
  struct Warning<S, true>
  {
    static inline void issue()
    {
      Dune::XT::Common::TimedLogger().get("gdt.l2localprojectionlocalizableoperator").warn()
          << "You are using this operator to project onto a continuous discrete function space (see below)!\n"
          << "Consider to use L2GlobalProjectionLocalizableOperator instead!\n"
          << "You can disable this warning by defining "
          << "DUNE_GDT_PROJECTIONS_L2_LOCAL_LOCALIZABLE_DISABLE_WARNING\n"
          << "at compile time or by disabling the Dune::XT::Common::TimedLogger() instance at runtime.\n"
          << "The type of the range space is: " << Dune::XT::Common::Typename<S>::value() << std::endl;
    } // ... issue_warning(...)
  };

  template <class S>
  static inline void issue_warning(const S&)
  {
#ifndef DUNE_GDT_PROJECTIONS_L2_LOCAL_LOCALIZABLE_DISABLE_WARNING
    Warning<S, S::continuous>::issue();
#endif
  }

  const LocalOperatorType local_operator_;
}; // class L2LocalProjectionLocalizableOperator


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Grid::is_layer<GridLayerType>::value && XT::Functions::is_localizable_function<SourceType>::value
                  && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2LocalProjectionLocalizableOperator<GridLayerType,
                                                                   SourceType,
                                                                   DiscreteFunction<SpaceType, VectorType>>>>::type
    make_local_l2_projection_localizable_operator(const GridLayerType& grid_layer,
                                                  const SourceType& source,
                                                  DiscreteFunction<SpaceType, VectorType>& range,
                                                  const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2LocalProjectionLocalizableOperator<GridLayerType,
                                                                            SourceType,
                                                                            DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, grid_layer, source, range);
} // ... make_local_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<L2LocalProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                   SourceType,
                                                                   DiscreteFunction<SpaceType, VectorType>>>>::type
    make_local_l2_projection_localizable_operator(const SourceType& source,
                                                  DiscreteFunction<SpaceType, VectorType>& range,
                                                  const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2LocalProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                            SourceType,
                                                                            DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, range.space().grid_layer(), source, range);
} // ... make_local_l2_projection_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2LocalProjectionOperator
    : public OperatorInterface<internal::L2LocalProjectionOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2LocalProjectionOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2LocalProjectionOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2LocalProjectionOperator(const size_t over_integrate, GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  L2LocalProjectionOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
  {
  }

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range) const
  {
    typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC> SourceType;
    L2LocalProjectionLocalizableOperator<GridLayerType, SourceType, DiscreteFunction<S, V>> op(
        over_integrate_, grid_layer_, source, range);
    op.apply();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template <class RangeType, class SourceType>
  void
  apply_inverse(const RangeType& /*range*/, SourceType& /*source*/, const XT::Common::Configuration& /*opts*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  std::vector<std::string> invert_options() const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

private:
  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class L2LocalProjectionOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2LocalProjectionOperator<GridLayerType>>>::type
make_local_l2_projection_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2LocalProjectionOperator<GridLayerType>>(over_integrate, grid_layer);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_L2_LOCAL_HH
