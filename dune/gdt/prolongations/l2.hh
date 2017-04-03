// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)

#ifndef DUNE_GDT_PROLONGATIONS_L2_HH
#define DUNE_GDT_PROLONGATIONS_L2_HH

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/discretefunction/default.hh>

#include "l2-local.hh"
#include "l2-global.hh"

namespace Dune {
namespace GDT {

// forward
template <class GridLayerImp, class FieldImp = double>
class L2ProlongationOperator;


namespace internal {


template <class GridLayerType, class SourceType, class RangeType>
class L2ProlongationLocalizableOperatorTraits
{
  static_assert(is_const_discrete_function<SourceType>::value, "");
  static_assert(is_discrete_function<RangeType>::value, "");

  template <class G, class S, class R, bool c = true>
  struct Helper
  {
    typedef L2GlobalProlongationLocalizableOperator<G, S, R> type;
  };

  template <class G, class S, class R>
  struct Helper<G, S, R, false>
  {
    typedef L2LocalProlongationLocalizableOperator<G, S, R> type;
  };

public:
  typedef typename Helper<GridLayerType, SourceType, RangeType, RangeType::SpaceType::continuous>::type BaseType;
}; // class L2ProlongationLocalizableOperatorTraits


template <class GridLayerType, class FieldImp = double>
class L2ProlongationOperatorTraits
{
public:
  typedef L2ProlongationOperator<GridLayerType, FieldImp> derived_type;
  typedef FieldImp FieldType;
};


} // namespace internal


template <class GridLayerImp, class SourceImp, class RangeImp>
class L2ProlongationLocalizableOperator
    : public internal::L2ProlongationLocalizableOperatorTraits<GridLayerImp, SourceImp, RangeImp>::BaseType
{
  typedef
      typename internal::L2ProlongationLocalizableOperatorTraits<GridLayerImp, SourceImp, RangeImp>::BaseType BaseType;

public:
  template <class... Args>
  explicit L2ProlongationLocalizableOperator(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


template <class GridLayerType, class SS, class SV, class RS, class RV>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2ProlongationLocalizableOperator<GridLayerType,
                                                                          ConstDiscreteFunction<SS, SV>,
                                                                          DiscreteFunction<RS, RV>>>>::type
make_l2_prolongation_localizable_operator(const GridLayerType& grid_layer,
                                          const ConstDiscreteFunction<SS, SV>& source,
                                          DiscreteFunction<RS, RV>& range,
                                          const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProlongationLocalizableOperator<GridLayerType,
                                                                         ConstDiscreteFunction<SS, SV>,
                                                                         DiscreteFunction<RS, RV>>>(
      over_integrate, grid_layer, source, range);
} // ... make_l2_prolongation_localizable_operator(...)

template <class SS, class SV, class RS, class RV>
std::unique_ptr<L2ProlongationLocalizableOperator<typename RS::GridLayerType,
                                                  ConstDiscreteFunction<SS, SV>,
                                                  DiscreteFunction<RS, RV>>>
make_l2_prolongation_localizable_operator(const ConstDiscreteFunction<SS, SV>& source,
                                          DiscreteFunction<RS, RV>& range,
                                          const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProlongationLocalizableOperator<typename RS::GridLayerType,
                                                                         ConstDiscreteFunction<SS, SV>,
                                                                         DiscreteFunction<RS, RV>>>(
      over_integrate, range.space().grid_layer(), source, range);
} // ... make_l2_prolongation_localizable_operator(...)


template <class GridLayerImp, class FieldImp>
class L2ProlongationOperator : public OperatorInterface<internal::L2ProlongationOperatorTraits<GridLayerImp, FieldImp>>
{
  typedef OperatorInterface<internal::L2ProlongationOperatorTraits<GridLayerImp, FieldImp>> BaseType;

public:
  typedef internal::L2ProlongationOperatorTraits<GridLayerImp, FieldImp> Traits;
  typedef GridLayerImp GridLayerType;
  using typename BaseType::FieldType;

private:
  using E = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype D;
  static const size_t d = GridLayerType::dimension;

public:
  L2ProlongationOperator(const size_t over_integrate, GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(over_integrate)
  {
  }

  L2ProlongationOperator(GridLayerType grid_layer)
    : grid_layer_(grid_layer)
    , over_integrate_(0)
  {
  }

  template <class SS, class SV, class RS, class RV>
  void apply(const ConstDiscreteFunction<SS, SV>& source, DiscreteFunction<RS, RV>& range) const
  {
    redirect<RS::continuous>::apply(grid_layer_, source, range, over_integrate_);
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
    static void apply(const GridLayerType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2GlobalProlongationLocalizableOperator<GridLayerType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  template <bool anything>
  struct redirect<false, anything>
  {
    template <class SourceType, class RangeType>
    static void apply(const GridLayerType& grd_vw, const SourceType& src, RangeType& rng, const size_t over_integrate)
    {
      L2LocalProlongationLocalizableOperator<GridLayerType, SourceType, RangeType>(over_integrate, grd_vw, src, rng)
          .apply();
    }
  };

  GridLayerType grid_layer_;
  const size_t over_integrate_;
}; // class L2ProlongationOperator


template <class GridLayerType>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value,
                        std::unique_ptr<L2ProlongationOperator<GridLayerType>>>::type
make_l2_prolongation_operator(const GridLayerType& grid_layer, const size_t over_integrate = 0)
{
  return Dune::XT::Common::make_unique<L2ProlongationOperator<GridLayerType>>(over_integrate, grid_layer);
}


template <class GridLayerType, class SS, class SV, class RS, class RV>
typename std::enable_if<XT::Grid::is_layer<GridLayerType>::value, void>::type
prolong_l2(const GridLayerType& grid_layer,
           const ConstDiscreteFunction<SS, SV>& source,
           DiscreteFunction<RS, RV>& range,
           const size_t over_integrate = 0)
{
  make_l2_prolongation_operator(grid_layer, over_integrate)->apply(source, range);
}

template <class SS, class SV, class RS, class RV>
void prolong_l2(const ConstDiscreteFunction<SS, SV>& source,
                DiscreteFunction<RS, RV>& range,
                const size_t over_integrate = 0)
{
  make_l2_prolongation_operator(range.space().grid_layer(), over_integrate)->apply(source, range);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_L2_HH
