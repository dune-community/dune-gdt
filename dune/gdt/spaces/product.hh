// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_PRODUCT_HH
#define DUNE_GDT_SPACES_PRODUCT_HH

#include <tuple>

#include <dune/stuff/common/tuple.hh>

#include <dune/gdt/basefunctionset/default/product.hh>
#include <dune/gdt/mapper/default/product.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {


template <class... SpaceImps>
class ProductSpace;


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
class ProductSpaceTraits
{
public:
  typedef ProductSpace<SpaceImps...> derived_type;
  typedef typename std::tuple_element<0, std::tuple<SpaceImps...>>::type::GridViewType GridViewType;
  static const size_t dimDomain    = GridViewType::dimension;
  static const size_t dimRange     = GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange;
  static const size_t dimRangeCols = 1;
  typedef typename GridViewType::IndexSet BackendType;
  typedef typename std::tuple_element<0, std::tuple<SpaceImps...>>::type::RangeFieldType RangeFieldType;
  typedef typename std::tuple<SpaceImps...> SpaceTupleType;
  typedef typename Dune::GDT::Mapper::DefaultProductMapper<GridViewType, typename SpaceImps::MapperType...> MapperType;
  typedef typename Dune::GDT::BaseFunctionSet::ProductBaseFunctionSet<typename SpaceImps::BaseFunctionSetType...>
      BaseFunctionSetType;
  static const int polOrder    = internal::maxPolOrder<SpaceImps...>::polOrder;
  static const bool continuous = internal::allContinuous<SpaceImps...>::value;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view                       = true;
  typedef CommunicationChooser<GridViewType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;
}; // class ProductSpaceTraits


} // namespace internal


template <class... SpaceImps>
class ProductSpace
    : public Dune::GDT::ProductSpaceInterface<internal::ProductSpaceTraits<SpaceImps...>,
                                              std::tuple_element<0, std::tuple<SpaceImps...>>::type::dimDomain,
                                              GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange, 1>
{
  typedef ProductSpace<SpaceImps...> ThisType;
  typedef Dune::GDT::ProductSpaceInterface<internal::ProductSpaceTraits<SpaceImps...>,
                                           std::tuple_element<0, std::tuple<SpaceImps...>>::type::dimDomain,
                                           GDT::BaseFunctionSet::internal::SumDimRange<SpaceImps...>::dimRange,
                                           1> BaseType;

public:
  using typename BaseType::Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;

public:
  using typename BaseType::CommunicatorType;
  using typename BaseType::SpaceTupleType;

  ProductSpace(SpaceImps&&... spaces)
    : spaces_(std::make_tuple(spaces...))
    , product_mapper_(spaces...)
    , communicator_(CommunicationChooserType::create(std::get<0>(spaces_).grid_view()))
  {
  }

  ProductSpace(const ThisType& other)
    : spaces_(other.spaces_)
    , product_mapper_(other.product_mapper_)
    , communicator_(CommunicationChooserType::create(std::get<0>(spaces_).grid_view()))
  {
  }

  ProductSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  // These methods are required by ProductSpaceInterface
  template <size_t ii>
  const typename std::tuple_element<ii, SpaceTupleType>::type& factor() const
  {
    return std::get<ii>(spaces_);
  }

  const MapperType& mapper() const
  {
    return product_mapper_;
  }

  // The remaining methods are redirected to Default
  const GridViewType& grid_view() const
  {
    return std::get<0>(spaces_).grid_view();
  }

  const BackendType& backend() const
  {
    return std::get<0>(spaces_).grid_view().indexSet();
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return base_function_set_helper(entity, typename DSC::create_indices<sizeof...(SpaceImps)>::type());
  }

  CommunicatorType& communicator() const
  {
    return *communicator_;
  }

private:
  template <size_t... S>
  BaseFunctionSetType base_function_set_helper(const EntityType& entity, DSC::indices<S...>) const
  {
    return BaseFunctionSetType(entity, std::get<S>(spaces_).base_function_set(entity)...);
  }

  std::tuple<SpaceImps...> spaces_;
  const MapperType product_mapper_;
  const std::unique_ptr<CommunicatorType> communicator_;
}; // class DefaultProductSpace


#include <dune/gdt/spaces/cg.hh>

template <class GridViewImp, int polynomialOrder, size_t domainDim, class RangeFieldImp>
class TaylorHoodSpace : public ProductSpace<CG::PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1>,
                                            CG::PdelabBased<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1>>
{
  typedef CG::PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, domainDim, 1> VelocitySpaceType;
  typedef CG::PdelabBased<GridViewImp, polynomialOrder - 1, RangeFieldImp, 1, 1> PressureSpaceType;
  typedef ProductSpace<VelocitySpaceType, PressureSpaceType> BaseType;

public:
  TaylorHoodSpace(const GridViewImp& grid_view)
    : BaseType(VelocitySpaceType(grid_view), PressureSpaceType(grid_view))
  {
  }
};


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PRODUCT_HH
