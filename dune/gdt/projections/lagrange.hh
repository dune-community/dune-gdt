// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROJECTIONS_LAGRANGE_HH
#define DUNE_GDT_PROJECTIONS_LAGRANGE_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/lagrange-projection.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \todo Add a check if the range space provides Lagrange points (after implementing the appropriate interface/mixin
 *       for those spaces).
 * \note Do we have to set all range DoFs to infinity here?
 */
template <class GridLayerImp, class SourceImp, class RangeImp>
class LagrangeProjectionLocalizableOperator : public LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp>
{
  typedef LocalizableOperatorBase<GridLayerImp, SourceImp, RangeImp> BaseType;

public:
  LagrangeProjectionLocalizableOperator(GridLayerImp grd_vw,
                                        const SourceImp& src,
                                        RangeImp& rng,
                                        const XT::Common::Parameter& param = {})
    : BaseType(grd_vw, src, rng)
    , local_operator_(param)
  {
    this->append(local_operator_);
  }

private:
  const LocalLagrangeProjectionOperator local_operator_;
}; // class LagrangeProjectionLocalizableOperator


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Grid::is_layer<GridLayerType>::value && XT::Functions::is_localizable_function<SourceType>::value
                  && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<LagrangeProjectionLocalizableOperator<GridLayerType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_lagrange_projection_localizable_operator(const GridLayerType& grid_layer,
                                                  const SourceType& source,
                                                  DiscreteFunction<SpaceType, VectorType>& range,
                                                  const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<LagrangeProjectionLocalizableOperator<GridLayerType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      grid_layer, source, range, param);
}


template <class SourceType, class SpaceType, class VectorType>
typename std::
    enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                  && XT::LA::is_vector<VectorType>::value,
              std::unique_ptr<LagrangeProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                    SourceType,
                                                                    DiscreteFunction<SpaceType, VectorType>>>>::type
    make_lagrange_projection_localizable_operator(const SourceType& source,
                                                  DiscreteFunction<SpaceType, VectorType>& range,
                                                  const XT::Common::Parameter& param = {})
{
  return Dune::XT::Common::make_unique<LagrangeProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                             SourceType,
                                                                             DiscreteFunction<SpaceType, VectorType>>>(
      range.space().grid_layer(), source, range, param);
}


// forward
template <class GridLayerImp, class FieldImp = double>
class LagrangeProjectionOperator;


namespace internal {


template <class GridLayerImp, class FieldImp>
class LagrangeProjectionOperatorTraits
{
public:
  typedef LagrangeProjectionOperator<GridLayerImp, FieldImp> derived_type;
  typedef NoJacobian JacobianType;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class FieldImp>
class LagrangeProjectionOperator
    : public OperatorInterface<internal::LagrangeProjectionOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::LagrangeProjectionOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::LagrangeProjectionOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  LagrangeProjectionOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
  {
  }

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range,
             const XT::Common::Parameter& param = {}) const
  {
    make_lagrange_projection_localizable_operator(grid_layer_, source, range, param)->apply();
  }

  template <class RangeType, class SourceType>
  FieldType
  apply2(const RangeType& /*range*/, const SourceType& /*source*/, const XT::Common::Parameter& param = {}) const
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
}; // class LagrangeProjectionOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<LagrangeProjectionOperator<GridLayerType>>>::type
make_lagrange_projection_operator(const GridLayerType& grid_layer)
{
  return Dune::XT::Common::make_unique<LagrangeProjectionOperator<GridLayerType>>(grid_layer);
}


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value
                            && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_lagrange(const GridLayerType& grid_layer,
                 const SourceType& source,
                 DiscreteFunction<SpaceType, VectorType>& range,
                 const XT::Common::Parameter& param = {})
{
  make_lagrange_projection_operator(grid_layer)->apply(source, range, param);
}


template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_lagrange(const SourceType& source,
                 DiscreteFunction<SpaceType, VectorType>& range,
                 const XT::Common::Parameter& param = {})
{
  make_lagrange_projection_operator(range.space().grid_layer())->apply(source, range, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_LAGRANGE_HH
