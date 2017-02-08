// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_FV_FvProductSpace_HH
#define DUNE_GDT_SPACES_FV_FvProductSpace_HH

#include <dune/xt/common/tuple.hh>

#include "default.hh"
#include "interface.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class FvProductSpace
{
  static_assert(Dune::AlwaysFalse<GridViewImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GridViewImp, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FvProductSpaceTraits : public FvSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef FvSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  typedef FvProductSpace<GridViewImp, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  using typename BaseType::GridViewType;
  static const size_t dimDomain = GridViewType::dimension;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  using typename BaseType::RangeFieldType;
  typedef typename Dune::GDT::FvSpace<GridViewType, RangeFieldType, 1, dimRangeCols> FactorSpaceType;
  typedef typename Dune::XT::Common::make_identical_tuple<FactorSpaceType, dimRange>::type SpaceTupleType;
  typedef typename Dune::GDT::FvProductMapper<GridViewType, dimRange, 1> MapperType;
}; // class FvProductSpaceTraits


} // namespace internal


template <class GridViewImp, class RangeFieldImp, size_t rangeDim>
class FvProductSpace<GridViewImp, RangeFieldImp, rangeDim, 1>
    : public Dune::GDT::FvSpaceInterface<internal::FvProductSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1>,
                                         GridViewImp::dimension,
                                         rangeDim,
                                         1>,
      public Dune::GDT::ProductSpaceInterface<internal::FvProductSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1>,
                                              GridViewImp::dimension,
                                              rangeDim,
                                              1>
{
  typedef FvProductSpace<GridViewImp, RangeFieldImp, rangeDim, 1> ThisType;
  typedef Dune::GDT::FvSpaceInterface<internal::FvProductSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1>,
                                      GridViewImp::dimension,
                                      rangeDim,
                                      1>
      BaseType;
  typedef FvSpace<GridViewImp, RangeFieldImp, rangeDim, 1> FvSpaceFVSpaceType;

public:
  typedef typename internal::FvProductSpaceTraits<GridViewImp, RangeFieldImp, rangeDim, 1> Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;
  typedef typename Traits::FactorSpaceType FactorSpaceType;

  FvProductSpace(GridViewType gv)
    : default_fv_space_(gv)
    , product_fv_mapper_(gv)
    , factor_space_(gv)
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
  const FvSpaceFVSpaceType default_fv_space_;
  const MapperType product_fv_mapper_;
  const FactorSpaceType factor_space_;
}; // class FvProductSpace< ..., 1 >


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_FV_FvProductSpace_HH
