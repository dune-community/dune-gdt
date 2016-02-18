// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_PLAYGROUND_SPACES_DG_PDELABPRODUCT_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DG_PDELABPRODUCT_HH

#include <tuple>

#include <dune/gdt/playground/mapper/productdgpdelab.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/dg/interface.hh>
#include <dune/gdt/playground/spaces/dg/pdelab.hh>

namespace Dune {
namespace GDT {
namespace Spaces {
namespace DG {


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class PdelabBasedProduct
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



template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols >
class PdelabBasedProductTraits
    : public PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
public:
  typedef PdelabBasedProduct< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  using typename BaseType::GridViewType;
  static const int polOrder = BaseType::polOrder;
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeFieldType;
  typedef Mapper::ProductDG< BackendType, rangeDim, rangeDimCols >          MapperType;
  using BaseType::part_view_type;
  using BaseType::needs_grid_view;

  typedef typename Dune::GDT::Spaces::DG::PdelabBased< GridViewType, polOrder, RangeFieldType, 1, rangeDimCols >  FactorSpaceType;
  typedef typename make_identical_tuple< FactorSpaceType, rangeDim >::type                                        SpaceTupleType;
};


} // namespace internal


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim >
class PdelabBasedProduct< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, 1 >
  : public Dune::GDT::ProductSpaceInterface< internal::PdelabBasedProductTraits< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, 1 >, GridViewImp::dimension, rangeDim, 1 >
  , public PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef PdelabBasedProduct< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, 1 >                          ThisType;
  typedef typename Dune::GDT::ProductSpaceInterface< internal::PdelabBasedProductTraits< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, 1 >, GridViewImp::dimension, rangeDim, 1 >  InterfaceType;
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > BaseType;
public:
  using typename InterfaceType::Traits;
  using typename InterfaceType::GridViewType;
  using typename InterfaceType::EntityType;
  using typename InterfaceType::BaseFunctionSetType;
  using typename InterfaceType::MapperType;
  using typename InterfaceType::CommunicatorType;
  using typename InterfaceType::BackendType;
  using InterfaceType::dimDomain;
  using InterfaceType::dimRange;
  using InterfaceType::dimRangeCols;
  typedef typename Traits::SpaceTupleType    SpaceTupleType;
  typedef typename Traits::FactorSpaceType   FactorSpaceType;

  PdelabBasedProduct(GridViewType gv)
    : BaseType(gv)
    , factor_space_(gv)
  {}

  PdelabBasedProduct(const ThisType& other)
    : BaseType(other)
    , factor_space_(other.factor_space_)
  {}

  PdelabBasedProduct(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  template< size_t ii >
  const FactorSpaceType& factor() const
  {
    return factor_space_;
  }

private:
    const FactorSpaceType factor_space_;
}; // class DefaultProduct< ..., 1 >


} // namespace DG
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_DG_PDELABPRODUCT_HH
