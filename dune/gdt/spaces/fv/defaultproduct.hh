// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH
#define DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH

#include <dune/stuff/common/tuple.hh>

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


template< class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols >
class DefaultProductTraits
    : public DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols >
{
  typedef DefaultTraits< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols >  BaseType;
public:
  typedef DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  using typename BaseType::GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  using typename BaseType::RangeFieldType;
  typedef typename Dune::GDT::Spaces::FV::Default< GridViewType, RangeFieldType, 1, dimRangeCols >  FactorSpaceType;
  typedef typename DSC::make_identical_tuple< FactorSpaceType, dimRange >::type                     SpaceTupleType;
  typedef typename Dune::GDT::Mapper::ProductFiniteVolume< GridViewType, dimRange, 1 >              MapperType;
}; // class DefaultProductTraits


} // namespace internal


template< class GridViewImp, class RangeFieldImp, size_t rangeDim >
class DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, 1 >
  : public Dune::GDT::Spaces::ProductFVInterface< internal::DefaultProductTraits< GridViewImp, RangeFieldImp, rangeDim, 1 >,
                                                  GridViewImp::dimension, rangeDim, 1 >
{
  typedef DefaultProduct< GridViewImp, RangeFieldImp, rangeDim, 1 >                  ThisType;
  typedef Dune::GDT::Spaces::ProductFVInterface
        < internal::DefaultProductTraits< GridViewImp, RangeFieldImp, rangeDim, 1 >,
          GridViewImp::dimension, rangeDim, 1 >                                      BaseType;
  typedef Default< GridViewImp, RangeFieldImp, rangeDim, 1 > DefaultFVSpaceType;
public:
  using typename BaseType::Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;
  typedef typename Traits::FactorSpaceType FactorSpaceType;

  DefaultProduct(GridViewType gv)
    : default_fv_space_(gv)
    , product_fv_mapper_(gv)
    , factor_space_(gv)
  {}

  DefaultProduct(const ThisType& other) = default;
  DefaultProduct(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  // These methods are required by ProductSpaceInterface
  template< size_t ii >
  const FactorSpaceType factor() const
  {
    return factor_space_;
  }

  const MapperType& mapper() const
  {
    return product_fv_mapper_;
  }

  // The remaining methods are redirected to Default
  const GridViewType& grid_view() const
  {
    return default_fv_space_.grid_view();
  }

  const BackendType& backend() const
  {
    return default_fv_space_.backend();
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return default_fv_space_.base_function_set(entity);
  }

  CommunicatorType& communicator() const
  {
    return default_fv_space_.communicator();
  }

private:
  const DefaultFVSpaceType default_fv_space_;
  const MapperType product_fv_mapper_;
  const FactorSpaceType factor_space_;
}; // class DefaultProduct< ..., 1 >


} // namespace FV
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_DEFAULTPRODUCT_HH
