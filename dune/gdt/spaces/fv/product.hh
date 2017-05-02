// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_FV_FvProductSpace_HH
#define DUNE_GDT_SPACES_FV_FvProductSpace_HH

#include <dune/xt/common/tuple.hh>

#include "default.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits and to allow for specialization
template <class GridLayerImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FvProductSpace
{
  static_assert(Dune::AlwaysFalse<GridLayerImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GridLayerImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FvProductSpaceTraits : public FvSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef FvSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  typedef FvProductSpace<GridLayerImp, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  using typename BaseType::GridLayerType;
  static const size_t dimDomain = GridLayerType::dimension;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  using typename BaseType::RangeFieldType;
  typedef typename Dune::GDT::FvSpace<GridLayerType, RangeFieldType, 1, dimRangeCols> FactorSpaceType;
  typedef typename Dune::XT::Common::make_identical_tuple<FactorSpaceType, dimRange>::type SpaceTupleType;
  typedef typename Dune::GDT::FvProductMapper<GridLayerType, dimRange, 1> MapperType;
}; // class FvProductSpaceTraits


} // namespace internal


template <class GridLayerImp, class RangeFieldImp, size_t rangeDim>
class FvProductSpace<GridLayerImp, RangeFieldImp, rangeDim, 1>
    : public Dune::GDT::FvSpaceInterface<internal::FvProductSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1>,
                                         GridLayerImp::dimension,
                                         rangeDim,
                                         1>,
      public Dune::GDT::ProductSpaceInterface<internal::FvProductSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1>,
                                              GridLayerImp::dimension,
                                              rangeDim,
                                              1>
{
  typedef FvProductSpace<GridLayerImp, RangeFieldImp, rangeDim, 1> ThisType;
  typedef Dune::GDT::FvSpaceInterface<internal::FvProductSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1>,
                                      GridLayerImp::dimension,
                                      rangeDim,
                                      1>
      BaseType;
  typedef FvSpace<GridLayerImp, RangeFieldImp, rangeDim, 1> FvSpaceFVSpaceType;

public:
  typedef typename internal::FvProductSpaceTraits<GridLayerImp, RangeFieldImp, rangeDim, 1> Traits;
  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;
  typedef typename Traits::FactorSpaceType FactorSpaceType;

  FvProductSpace(GridLayerType grd_layr)
    : default_fv_space_(grd_layr)
    , product_fv_mapper_(grd_layr)
    , factor_space_(grd_layr)
  {
  }

  FvProductSpace(const ThisType& other) = default;
  FvProductSpace(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  // These methods are required by ProductSpaceInterface
  template <size_t ii>
  const FactorSpaceType& factor() const
  {
    return factor_space_;
  }

  const MapperType& mapper() const
  {
    return product_fv_mapper_;
  }

  // The remaining methods are redirected to FvSpace
  const GridLayerType& grid_layer() const
  {
    return default_fv_space_.grid_layer();
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
  const FvSpaceFVSpaceType default_fv_space_;
  const MapperType product_fv_mapper_;
  const FactorSpaceType factor_space_;
}; // class FvProductSpace< ..., 1 >


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_FvProductSpace_HH
