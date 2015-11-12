// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH
#define DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH

#include <tuple>

#include <dune/gdt/mapper/default/productfv.hh>

#include "default.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace FV {


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class DefaultProduct
{
  static_assert(Dune::AlwaysFalse< GridViewImp >::value, "Untested for these dimensions!");
};


namespace internal {

// from https://stackoverflow.com/questions/16853552/how-to-create-a-type-list-for-variadic-templates-that-contains-n-times-the-sam

// in the end, we would like to have something like indices< 1, 2, 3 > for N = 3
template< std::size_t... >
struct indices {};

// we want to call this with empty Indices, i.e. create_indices< N >::type == indices< 1, 2, 3 > for N = 3
template< std::size_t N, std::size_t... Indices>
struct create_indices : create_indices< N-1, N-1, Indices...> {};

// terminating template
template< std::size_t... Indices >
struct create_indices< 0, Indices...> {
  typedef indices<Indices...> type;
};

// T_aliased< T, Index > is always the type T, no matter what Index is
template<typename T, std::size_t index>
using T_aliased = T;

// make_identical_tuple< T, N >::type is a std::tuple< T, ... , T > with a length of N
template< typename T, std::size_t N, typename I = typename create_indices< N >::type >
struct make_identical_tuple;

template< typename T, std::size_t N, std::size_t ...Indices >
struct make_identical_tuple< T, N, indices< Indices... > >
{
    using type = std::tuple<T_aliased<T, Indices>...>;
};



template< class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols >
class DefaultProductTraits
    : public DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols >  BaseType;
public:
  typedef DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  using typename BaseType::GridViewType;
  using typename BaseType::EntityType;
  static const size_t dimDomain = GridViewType::dimension;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  using typename BaseType::RangeFieldType;
  typedef typename Dune::GDT::Spaces::FV::Default< GridViewType, RangeFieldType, 1, dimRangeCols >  FactorSpaceType;
  typedef typename make_identical_tuple< FactorSpaceType, dimRange >::type                          SpaceTupleType;
  typedef typename Dune::GDT::Mapper::ProductFiniteVolume< GridViewType, dimRange, 1 >              FactorMapperType;
}; // class DefaultProductTraits


} // namespace internal


template< class GridViewImp, class RangeFieldImp, size_t rangeDim >
class DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, 1 >
  : public Dune::GDT::Spaces::ProductFVInterface< internal::DefaultProductTraits< GridViewImp, RangeFieldImp, rangeDim, 1 > >
{
  typedef DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, 1 >                          ThisType;
  typedef Dune::GDT::Spaces::ProductFVInterface
        < internal::DefaultProductTraits< GridViewImp, RangeFieldImp, rangeDim, 1 > >        BaseType;
public:
  using typename BaseType::Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::RangeFieldType;
private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
public:
  using typename BaseType::CommunicatorType;
  typedef typename Traits::FactorMapperType FactorMapperType;
  typedef typename Traits::SpaceTupleType   SpaceTupleType;
  typedef typename Traits::FactorSpaceType  FactorSpaceType;

  DefaultProduct(GridViewType gv)
    : grid_view_(gv)
    , mapper_(grid_view_)
    , communicator_(CommunicationChooserType::create(grid_view_))
    , factor_space_(grid_view_)
    , factor_mapper_(grid_view_)
  {}

  DefaultProduct(const ThisType& other)
    : grid_view_(other.grid_view_)
    , mapper_(other.mapper_)
    , communicator_(CommunicationChooserType::create(grid_view_))
    , factor_space_(other.factor_space_)
    , factor_mapper_(other.factor_mapper_)
  {}

  DefaultProduct(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const FactorMapperType& factor_mapper() const
  {
    return factor_mapper_;
  }

  template< size_t ii >
  const FactorSpaceType factor() const
  {
    return factor_space_;
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return grid_view_.indexSet();
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(entity);
  }

  CommunicatorType& communicator() const
  {
    // no need to prepare the communicator, since we are not pdelab based
    return *communicator_;
  }

private:
  const GridViewType grid_view_;
  const MapperType mapper_;
  const std::unique_ptr< CommunicatorType > communicator_;
  const FactorSpaceType factor_space_;
  const FactorMapperType factor_mapper_;
}; // class DefaultProduct< ..., 1 >


} // namespace FV
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH
