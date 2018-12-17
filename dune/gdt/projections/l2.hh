// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_PROJECTIONS_L2_HH
#define DUNE_GDT_PROJECTIONS_L2_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/type_traits.hh>

#include "l2-local.hh"
#include "l2-global.hh"

namespace Dune {
namespace GDT {

// forward
template <class GridLayerImp, class FieldImp = double>
class L2ProjectionOperator;


namespace internal {


template <class GridLayerType, class SourceType, class RangeType>
class L2ProjectionLocalizableOperatorTraits
{
  static_assert(is_discrete_function<RangeType>::value, "");

  template <class G, class S, class R, bool c = true>
  struct Helper
  {
    typedef L2GlobalProjectionLocalizableOperator<G, S, R> type;
  };

  template <class G, class S, class R>
  struct Helper<G, S, R, false>
  {
    typedef L2LocalProjectionLocalizableOperator<G, S, R> type;
  };

public:
  typedef typename Helper<GridLayerType, SourceType, RangeType, RangeType::SpaceType::continuous>::type BaseType;
}; // class L2ProjectionLocalizableOperatorTraits


template <class GridLayerType, class FieldImp = double>
class L2ProjectionOperatorTraits
{
public:
  typedef L2ProjectionOperator<GridLayerType, FieldImp> derived_type;
  typedef FieldImp FieldType;
  typedef NoJacobian JacobianType;
};


} // namespace internal


template <class GridLayerImp, class SourceImp, class RangeImp>
class L2ProjectionLocalizableOperator
  : public internal::L2ProjectionLocalizableOperatorTraits<GridLayerImp, SourceImp, RangeImp>::BaseType
{
  typedef
      typename internal::L2ProjectionLocalizableOperatorTraits<GridLayerImp, SourceImp, RangeImp>::BaseType BaseType;

public:
  template <class... Args>
  explicit L2ProjectionLocalizableOperator(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {}
};


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<
    XT::Grid::is_layer<GridLayerType>::value && XT::Functions::is_localizable_function<SourceType>::value
        && is_space<SpaceType>::value && XT::LA::is_vector<VectorType>::value,
    std::unique_ptr<
        L2ProjectionLocalizableOperator<GridLayerType, SourceType, DiscreteFunction<SpaceType, VectorType>>>>::type
make_l2_projection_localizable_operator(const GridLayerType& grid_layer,
                                        const SourceType& source,
                                        DiscreteFunction<SpaceType, VectorType>& range,
                                        const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<
      L2ProjectionLocalizableOperator<GridLayerType, SourceType, DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, grid_layer, source, range);
} // ... make_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        std::unique_ptr<L2ProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                        SourceType,
                                                                        DiscreteFunction<SpaceType, VectorType>>>>::type
make_l2_projection_localizable_operator(const SourceType& source,
                                        DiscreteFunction<SpaceType, VectorType>& range,
                                        const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProjectionLocalizableOperator<typename SpaceType::GridLayerType,
                                                                       SourceType,
                                                                       DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, range.space().grid_layer(), source, range);
} // ... make_l2_projection_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2ProjectionOperator : public OperatorInterface<internal::L2ProjectionOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2ProjectionOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2ProjectionOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2ProjectionOperator(const size_t over_integrate, GridLayerType grid_layer, const bool use_tbb = false)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
    , use_tbb_(use_tbb)
  {}

  L2ProjectionOperator(GridLayerType grid_layer, const bool use_tbb = false)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
    , use_tbb_(use_tbb)
  {}

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range,
             const XT::Common::Parameter& param = {}) const
  {
    redirect<S::continuous>::apply(grid_layer_, source, range, over_integrate_, param, use_tbb_);
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
  template <bool continuous = true, bool anything = true>
  struct redirect
  {
    template <class SourceType, class RangeType>
    static void apply(const GridLayerType& grd_vw,
                      const SourceType& src,
                      RangeType& rng,
                      const size_t over_integrate,
                      const XT::Common::Parameter& param = {},
                      const bool /*use_tbb*/ = false)
    {
      L2GlobalProjectionLocalizableOperator<GridLayerType, SourceType, RangeType>(
          over_integrate, grd_vw, src, rng, param)
          .apply(/*use_tbb*/);
    }
  };

  template <bool anything>
  struct redirect<false, anything>
  {
    template <class SourceType, class RangeType>
    static void apply(const GridLayerType& grd_vw,
                      const SourceType& src,
                      RangeType& rng,
                      const size_t over_integrate,
                      const XT::Common::Parameter& param = {},
                      const bool use_tbb = false)
    {
      L2LocalProjectionLocalizableOperator<GridLayerType, SourceType, RangeType>(
          over_integrate, grd_vw, src, rng, param)
          .apply(use_tbb);
    }
  };

  GridLayerType grid_layer_;
  const size_t over_integrate_;
  const bool use_tbb_;
}; // class L2ProjectionOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2ProjectionOperator<GridLayerType>>>::type
make_l2_projection_operator(const GridLayerType& grid_layer,
                            const size_t over_integrate = 0,
                            const bool use_tbb = false)
{
  return Dune::XT::Common::make_unique<L2ProjectionOperator<GridLayerType>>(over_integrate, grid_layer, use_tbb);
}


template <class GridLayerType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_l2(const GridLayerType& grid_layer,
           const SourceType& source,
           DiscreteFunction<SpaceType, VectorType>& range,
           const size_t over_integrate = 0,
           const XT::Common::Parameter& param = {},
           const bool use_tbb = false)
{
  make_l2_projection_operator(grid_layer, over_integrate, use_tbb)->apply(source, range, param);
}


template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_l2(const SourceType& source,
           DiscreteFunction<SpaceType, VectorType>& range,
           const size_t over_integrate = 0,
           const XT::Common::Parameter& param = {},
           const bool use_tbb = false)
{
  make_l2_projection_operator(range.space().grid_layer(), over_integrate, use_tbb)->apply(source, range, param);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_L2_HH
