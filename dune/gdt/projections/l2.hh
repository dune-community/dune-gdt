// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_PROJECTIONS_L2_HH
#define DUNE_GDT_PROJECTIONS_L2_HH

#include <dune/gdt/discretefunction/default.hh>

#include "l2-local.hh"
#include "l2-global.hh"

namespace Dune {
namespace GDT {

// forward
template <class GridViewImp, class FieldImp = double>
class L2ProjectionOperator;


namespace internal {


template <class GridViewType, class SourceType, class RangeType>
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
  typedef typename Helper<GridViewType, SourceType, RangeType, RangeType::SpaceType::continuous>::type BaseType;
}; // class L2ProjectionLocalizableOperatorTraits


template <class GridViewType, class FieldImp = double>
class L2ProjectionOperatorTraits
{
public:
  typedef L2ProjectionOperator<GridViewType, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridViewImp, class SourceImp, class RangeImp>
class L2ProjectionLocalizableOperator
    : public internal::L2ProjectionLocalizableOperatorTraits<GridViewImp, SourceImp, RangeImp>::BaseType
{
  typedef typename internal::L2ProjectionLocalizableOperatorTraits<GridViewImp, SourceImp, RangeImp>::BaseType BaseType;

public:
  template <class... Args>
  explicit L2ProjectionLocalizableOperator(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


template <class GridViewType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        std::unique_ptr<L2ProjectionLocalizableOperator<GridViewType, SourceType,
                                                                        DiscreteFunction<SpaceType, VectorType>>>>::type
make_l2_projection_localizable_operator(const GridViewType& grid_view, const SourceType& source,
                                        DiscreteFunction<SpaceType, VectorType>& range, const size_t over_integrate = 0)
{
  return Dune::XT::Common::
      make_unique<L2ProjectionLocalizableOperator<GridViewType, SourceType, DiscreteFunction<SpaceType, VectorType>>>(
          over_integrate, grid_view, source, range);
} // ... make_l2_projection_localizable_operator(...)

template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        std::unique_ptr<L2ProjectionLocalizableOperator<typename SpaceType::GridViewType, SourceType,
                                                                        DiscreteFunction<SpaceType, VectorType>>>>::type
make_l2_projection_localizable_operator(const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range,
                                        const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProjectionLocalizableOperator<typename SpaceType::GridViewType,
                                                                       SourceType,
                                                                       DiscreteFunction<SpaceType, VectorType>>>(
      over_integrate, range.space().grid_view(), source, range);
} // ... make_l2_projection_localizable_operator(...)


template <class GridViewImp, class FieldImp>
class L2ProjectionOperator : public OperatorInterface<internal::L2ProjectionOperatorTraits<GridViewImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2ProjectionOperatorTraits<GridViewImp, FieldImp>> BaseType;

public:
  typedef internal::L2ProjectionOperatorTraits<GridViewImp, FieldImp> Traits;
  typedef GridViewImp GridViewType;
  using typename BaseType::FieldType;

private:
  typedef typename XT::Grid::Entity<GridViewType>::Type E;
  typedef typename GridViewType::ctype D;
  static const size_t d = GridViewType::dimension;

public:
  L2ProjectionOperator(const size_t over_integrate, GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {
  }

  L2ProjectionOperator(GridViewType grid_view)
    : grid_view_(grid_view)
    , over_integrate_(0)
  {
  }

  template <class R, size_t r, size_t rC, class S, class V>
  void apply(const XT::Functions::LocalizableFunctionInterface<E, D, d, R, r, rC>& source,
             DiscreteFunction<S, V>& range) const
  {
    redirect<S::continuous>::apply(grid_view_, source, range, over_integrate_);
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& /*range*/, const SourceType& /*source*/) const
  {
    DUNE_THROW(NotImplemented, "Go ahead if you think this makes sense!");
  }

  template <class RangeType, class SourceType>
  void apply_inverse(const RangeType& /*range*/, SourceType& /*source*/,
                     const XT::Common::Configuration& /*opts*/) const
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
    static void apply(const GridViewType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2GlobalProjectionLocalizableOperator<GridViewType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  template <bool anything>
  struct redirect<false, anything>
  {
    template <class SourceType, class RangeType>
    static void apply(const GridViewType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2LocalProjectionLocalizableOperator<GridViewType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  GridViewType grid_view_;
  const size_t over_integrate_;
}; // class L2ProjectionOperator


template <class GridViewType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value,
                        std::unique_ptr<L2ProjectionOperator<GridViewType>>>::type
make_l2_projection_operator(const GridViewType& grid_view, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProjectionOperator<GridViewType>>(over_integrate, grid_view);
}


template <class GridViewType, class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Grid::is_layer<GridViewType>::value
                            && XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_l2(const GridViewType& grid_view, const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range,
           const size_t over_integrate = 0)
{
  make_l2_projection_operator(grid_view, over_integrate)->apply(source, range);
}


template <class SourceType, class SpaceType, class VectorType>
typename std::enable_if<XT::Functions::is_localizable_function<SourceType>::value && is_space<SpaceType>::value
                            && XT::LA::is_vector<VectorType>::value,
                        void>::type
project_l2(const SourceType& source, DiscreteFunction<SpaceType, VectorType>& range, const size_t over_integrate = 0)
{
  make_l2_projection_operator(range.space().grid_view(), over_integrate)->apply(source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROJECTIONS_L2_HH
