// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_SPACES_PRODUCT_HH
#define DUNE_GDT_SPACES_PRODUCT_HH

#include <tuple>

#include <dune/xt/common/tuple.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/basefunctionset/product.hh>
#include <dune/gdt/spaces/mapper/product.hh>
#include <dune/gdt/spaces/cg.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class... SpaceImps>
class DefaultProductSpace;


namespace internal {


template <class T>
constexpr T staticMax(T a, T b)
{
  return a > b ? a : b;
}

template <class FirstSpaceType, class... SpaceTypes>
struct maxPolOrder
{
  static const int polOrder = staticMax(FirstSpaceType::polOrder, maxPolOrder<SpaceTypes...>::polOrder);
};

template <class LastSpaceType>
struct maxPolOrder<LastSpaceType>
{
  static const int polOrder = LastSpaceType::polOrder;
};

template <class FirstSpaceType, class... SpaceTypes>
struct allContinuous
{
  static const bool value = FirstSpaceType::continuous ? allContinuous<SpaceTypes...>::value : false;
};

template <class LastSpaceType>
struct allContinuous<LastSpaceType>
{
  static const bool value = LastSpaceType::continuous;
};


template <class... SpaceImps>
class DefaultProductSpaceTraits
{
public:
  typedef DefaultProductSpace<SpaceImps...> derived_type;
  typedef typename XT::Common::tuple_element<0, std::tuple<SpaceImps...>>::type::GridLayerType GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  static const size_t dimRange = GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename GridLayerType::IndexSet BackendType;
  typedef typename XT::Common::tuple_element<0, std::tuple<SpaceImps...>>::type::RangeFieldType RangeFieldType;
  typedef typename std::tuple<SpaceImps...> SpaceTupleType;
  typedef typename Dune::GDT::DefaultProductMapper<GridLayerType, typename SpaceImps::MapperType...> MapperType;
  typedef typename Dune::GDT::BaseFunctionSet::ProductDefault<typename SpaceImps::BaseFunctionSetType...>
      BaseFunctionSetType;
  static const int polOrder = internal::maxPolOrder<SpaceImps...>::polOrder;
  static const bool continuous = internal::allContinuous<SpaceImps...>::value;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const bool needs_grid_view = true;
  typedef DofCommunicationChooser<GridLayerType> DofCommunicationChooserType;
  typedef typename DofCommunicationChooserType::Type DofCommunicatorType;
  static const constexpr Backends backend_type{std::tuple_element<0, std::tuple<SpaceImps...>>::Traits::backend_type};
}; // class ProductSpaceTraits


} // namespace internal


template <class... SpaceImps>
class DefaultProductSpace
  : public Dune::GDT::SpaceInterface<internal::DefaultProductSpaceTraits<SpaceImps...>,
                                     XT::Common::tuple_element<0, std::tuple<SpaceImps...>>::type::dimDomain,
                                     GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange,
                                     1>
  , public Dune::GDT::ProductSpaceInterface<internal::DefaultProductSpaceTraits<SpaceImps...>,
                                            XT::Common::tuple_element<0, std::tuple<SpaceImps...>>::type::dimDomain,
                                            GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange,
                                            1>
{
  typedef DefaultProductSpace<SpaceImps...> ThisType;
  typedef Dune::GDT::SpaceInterface<internal::DefaultProductSpaceTraits<SpaceImps...>,
                                    XT::Common::tuple_element<0, std::tuple<SpaceImps...>>::type::dimDomain,
                                    GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange,
                                    1>
      BaseType;

public:
  typedef typename internal::DefaultProductSpaceTraits<SpaceImps...> Traits;
  using typename BaseType::BackendType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::EntityType;
  using typename BaseType::GridLayerType;
  using typename BaseType::MapperType;

private:
  typedef typename Traits::DofCommunicationChooserType DofCommunicationChooserType;

public:
  using typename BaseType::DofCommunicatorType;
  typedef typename Traits::SpaceTupleType SpaceTupleType;

  DefaultProductSpace(const SpaceImps&... spaces)
    : spaces_(std::make_tuple(spaces...))
    , product_mapper_(spaces_)
    , communicator_(DofCommunicationChooserType::create(std::get<0>(spaces_).grid_layer()))
  {}

  DefaultProductSpace(const ThisType& other)
    : spaces_(other.spaces_)
    , product_mapper_(other.product_mapper_)
    , communicator_(DofCommunicationChooserType::create(std::get<0>(spaces_).grid_layer()))
  {}

  DefaultProductSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  // These methods are required by ProductSpaceInterface
  template <size_t ii>
  const typename XT::Common::tuple_element<ii, SpaceTupleType>::type& factor() const
  {
    return std::get<ii>(spaces_);
  }

  const MapperType& mapper() const
  {
    return product_mapper_;
  }

  // The remaining methods are redirected to Default
  const GridLayerType& grid_layer() const
  {
    return std::get<0>(spaces_).grid_layer();
  }

  GridLayerType& grid_layer()
  {
    return std::get<0>(spaces_).grid_layer();
  }

  const BackendType& backend() const
  {
    return std::get<0>(spaces_).grid_layer().indexSet();
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return base_function_set_helper(entity,
                                    typename Dune::XT::Common::make_index_sequence<sizeof...(SpaceImps)>::type());
  }

  DofCommunicatorType& dof_communicator() const
  {
    return *communicator_;
  }

  static constexpr bool associates_data_with(int codim)
  {
    return BaseType::associates_data_with(codim);
  }

private:
  template <size_t... S>
  BaseFunctionSetType base_function_set_helper(const EntityType& entity, Dune::XT::Common::index_sequence<S...>) const
  {
    return BaseFunctionSetType(entity, std::get<S>(spaces_).base_function_set(entity)...);
  }

  std::tuple<SpaceImps...> spaces_;
  const MapperType product_mapper_;
  const std::unique_ptr<DofCommunicatorType> communicator_;
}; // class DefaultProductSpace


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PRODUCT_HH
