// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH

#include <memory>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/type.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/constraints/conforming.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/xt/common/tuple.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/parallel.hh>
#include <dune/gdt/spaces/dg/interface.hh>
#include <dune/gdt/spaces/basefunctionset/dune-pdelab-wrapper.hh>
#include <dune/gdt/playground/spaces/mapper/dune-pdelab-wrapper.hh>

namespace Dune {
namespace GDT {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridLayerImp>::value, "Untested for this combination of dimensions!");
};


// forward, to be used in the traits and to allow for specialization
template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgProductSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridLayerImp>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapperTraits
{
public:
  typedef DunePdelabDgSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridLayerImp GridLayerType;
  static const int polOrder = polynomialOrder;
  static const bool continuous = false;

private:
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = GridLayerType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;

private:
  template <class G, bool single_geom, bool is_simplex, bool is_cube>
  struct FeMap
  {
    static_assert(Dune::AlwaysFalse<G>::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template <class G>
  struct FeMap<G, true, true, false>
  {
    static_assert(Dune::AlwaysFalse<G>::value, "Not yet implemented for simplicial grids!");
  };
  template <class G>
  struct FeMap<G, true, false, true>
  {
    typedef PDELab::QkDGLocalFiniteElementMap<DomainFieldType, RangeFieldType, polOrder, dimDomain> Type;
  };
  typedef typename GridLayerType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType<GridType>::v;
  static const bool simplicial_ =
      (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId == Impl::SimplexTopology<dimDomain>::type::id);
  static const bool cubic_ =
      (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId == Impl::CubeTopology<dimDomain>::type::id);
  typedef typename FeMap<GridType, single_geom_, simplicial_, cubic_>::Type FEMapType;

public:
  typedef PDELab::GridFunctionSpace<GridLayerType, FEMapType, PDELab::OverlappingConformingDirichletConstraints>
      BackendType;
  typedef DunePdelabDgMapperWrapper<BackendType> MapperType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef BaseFunctionSet::
      DunePdelabWrapper<BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, rangeDim, rangeDimCols>
          BaseFunctionSetType;
  static const XT::Grid::Backends part_view_type = XT::Grid::Backends::view;
  static const bool needs_grid_view = true;
  typedef CommunicationChooser<GridLayerType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;

private:
  friend class DunePdelabDgSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class DunePdelabDgSpaceWrapperTraits


template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class DunePdelabDgProductSpaceWrapperTraits
    : public DunePdelabDgSpaceWrapperTraits<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>
{
  typedef DunePdelabDgSpaceWrapperTraits<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> BaseType;

public:
  typedef DunePdelabDgProductSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>
      derived_type;
  using typename BaseType::GridLayerType;
  static const int polOrder = BaseType::polOrder;
  static const size_t dimDomain = GridLayerType::dimension;
  static const size_t dimRange = rangeDim;
  static const size_t dimRangeCols = rangeDimCols;
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeFieldType;
  typedef ProductDgMapper<BackendType, rangeDim, rangeDimCols> MapperType;
  using BaseType::part_view_type;
  using BaseType::needs_grid_view;

  typedef typename Dune::GDT::DunePdelabDgSpaceWrapper<GridLayerType, polOrder, RangeFieldType, 1, dimRangeCols>
      FactorSpaceType;
  typedef typename Dune::XT::Common::make_identical_tuple<FactorSpaceType, dimRange>::type SpaceTupleType;
};


} // namespace internal


template <class GridLayerImp, int polynomialOrder, class RangeFieldImp>
class DunePdelabDgSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public DgSpaceInterface<internal::
                                  DunePdelabDgSpaceWrapperTraits<GridLayerImp, polynomialOrder, RangeFieldImp, 1, 1>,
                              GridLayerImp::dimension,
                              1,
                              1>
{
  typedef DunePdelabDgSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;
  typedef DgSpaceInterface<internal::DunePdelabDgSpaceWrapperTraits<GridLayerImp, polynomialOrder, RangeFieldImp, 1, 1>,
                           GridLayerImp::dimension,
                           1,
                           1>
      BaseType;

public:
  using typename BaseType::Traits;

  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
  typedef typename Traits::FEMapType FEMapType;

public:
  using typename BaseType::CommunicatorType;

  DunePdelabDgSpaceWrapper(GridLayerType grd_vw)
    : grid_layer_(grd_vw)
    , fe_map_()
    , backend_(grid_layer_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridLayerImp>::create(grid_layer_))
    , communicator_prepared_(false)
  {
  }

  /**
   * \brief Copy ctor.
   * \note  Manually implemented bc of the std::mutex + communicator_  unique_ptr
   */
  DunePdelabDgSpaceWrapper(const ThisType& other)
    : grid_layer_(other.grid_layer_)
    , fe_map_()
    , backend_(grid_layer_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridLayerImp>::create(grid_layer_))
    , communicator_prepared_(false)
  {
    // make sure our new communicator is prepared if other's was
    if (other.communicator_prepared_)
      const auto& comm DUNE_UNUSED = this->communicator();
  }

  /**
   * \brief Move ctor.
   * \note  Manually implemented bc of the std::mutex.
   */
  DunePdelabDgSpaceWrapper(ThisType&& source)
    : grid_layer_(source.grid_layer_)
    , fe_map_(source.fe_map_)
    , backend_(source.backend_)
    , mapper_(source.mapper_)
    , communicator_(std::move(source.communicator_))
    , communicator_prepared_(source.communicator_prepared_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(backend_, entity);
  }

  CommunicatorType& communicator() const
  {
    DUNE_UNUSED std::lock_guard<std::mutex> gg(communicator_mutex_);
    if (!communicator_prepared_)
      communicator_prepared_ = CommunicationChooser<GridLayerType>::prepare(*this, *communicator_);
    return *communicator_;
  } // ... communicator(...)

private:
  GridLayerType grid_layer_;
  const FEMapType fe_map_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class DunePdelabDgSpaceWrapper< ..., 1 >


template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim>
class DunePdelabDgProductSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, 1>
    : public Dune::GDT::SpaceInterface<internal::DunePdelabDgProductSpaceWrapperTraits<GridLayerImp,
                                                                                       polynomialOrder,
                                                                                       RangeFieldImp,
                                                                                       rangeDim,
                                                                                       1>,
                                       GridLayerImp::dimension,
                                       rangeDim,
                                       1>,
      public Dune::GDT::ProductSpaceInterface<internal::DunePdelabDgProductSpaceWrapperTraits<GridLayerImp,
                                                                                              polynomialOrder,
                                                                                              RangeFieldImp,
                                                                                              rangeDim,
                                                                                              1>,
                                              GridLayerImp::dimension,
                                              rangeDim,
                                              1>
{
  typedef DunePdelabDgProductSpaceWrapper<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, 1> ThisType;
  typedef typename Dune::GDT::SpaceInterface<internal::DunePdelabDgProductSpaceWrapperTraits<GridLayerImp,
                                                                                             polynomialOrder,
                                                                                             RangeFieldImp,
                                                                                             rangeDim,
                                                                                             1>,
                                             GridLayerImp::dimension,
                                             rangeDim,
                                             1>
      BaseType;

public:
  typedef typename internal::
      DunePdelabDgProductSpaceWrapperTraits<GridLayerImp, polynomialOrder, RangeFieldImp, rangeDim, 1>
          Traits;
  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::MapperType;
  using typename BaseType::CommunicatorType;
  using typename BaseType::BackendType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  typedef typename Traits::SpaceTupleType SpaceTupleType;
  typedef typename Traits::FactorSpaceType FactorSpaceType;

  DunePdelabDgProductSpaceWrapper(GridLayerType grd_vw)
    : grid_layer_(grd_vw)
    , factor_space_(grid_layer_)
    , factor_mapper_(factor_space_.backend())
    , communicator_(CommunicationChooser<GridLayerImp>::create(grid_layer_))
    , communicator_prepared_(false)
  {
  }

  DunePdelabDgProductSpaceWrapper(const ThisType& other)
    : grid_layer_(other.grid_layer_)
    , factor_space_(other.factor_space_)
    , factor_mapper_(other.factor_mapper_)
    , communicator_(CommunicationChooser<GridLayerImp>::create(grid_layer_))
    , communicator_prepared_(false)
  {
    // make sure our new communicator is prepared if other's was
    if (other.communicator_prepared_)
      const auto& comm DUNE_UNUSED = this->communicator();
  }

  DunePdelabDgProductSpaceWrapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_layer_;
  }

  GridLayerType& grid_layer()
  {
    return grid_layer_;
  }

  const BackendType& backend() const
  {
    return factor_space_.backend();
  }

  const MapperType& mapper() const
  {
    return factor_mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(backend(), entity);
  }

  CommunicatorType& communicator() const
  {
    DUNE_UNUSED std::lock_guard<std::mutex> gg(communicator_mutex_);
    if (!communicator_prepared_)
      communicator_prepared_ = CommunicationChooser<GridLayerType>::prepare(*this, *communicator_);
    return *communicator_;
  } // ... communicator(...)

  template <size_t ii>
  const FactorSpaceType& factor() const
  {
    return factor_space_;
  }

private:
  const GridLayerType grid_layer_;
  const FactorSpaceType factor_space_;
  const MapperType factor_mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class DefaultProduct< ..., 1 >


#else // HAVE_DUNE_PDELAB


template <class GridLayerImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridLayerImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH
