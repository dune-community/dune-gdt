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

#ifndef DUNE_GDT_SPACES_RT_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_SPACES_RT_DUNE_PDELAB_WRAPPER_HH

#include <type_traits>
#include <limits>
#include <mutex>

#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>

#include <dune/gdt/spaces/basefunctionset/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/mapper/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabRtSpaceWrapper
{
  static_assert(AlwaysFalse<GridViewImp>::value, "Untested for these dimensions or polynomial order!");
}; // class DunePdelabRtSpaceWrapper


namespace internal {


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class DunePdelabRtSpaceWrapperTraits
{
  static_assert(XT::Grid::is_view<GridViewImp>::value, "");

public:
  typedef DunePdelabRtSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridViewImp GridLayerType;
  static const int polOrder = polynomialOrder;
  static const bool continuous = false;
  static_assert(polOrder == 0, "Untested!");
  static_assert(rangeDim == GridLayerType::dimension, "Untested!");
  static_assert(rangeDimCols == 1, "Untested!");
  static const constexpr Backends backend_type{Backends::pdelab};

private:
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = GridLayerType::dimension;

public:
  typedef RangeFieldImp RangeFieldType;

private:
  template <class G, bool single_geom, bool is_simplex, bool is_cube>
  struct FeMap
  {
    static_assert(AlwaysFalse<G>::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template <class G>
  struct FeMap<G, true, true, false>
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap<GridLayerType,
                                                       DomainFieldType,
                                                       RangeFieldType,
                                                       polOrder,
                                                       Dune::GeometryType::simplex>
        Type;
  };
  template <class G>
  struct FeMap<G, true, false, true>
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap<GridLayerType,
                                                       DomainFieldType,
                                                       RangeFieldType,
                                                       polOrder,
                                                       Dune::GeometryType::cube>
        Type;
  };
  typedef XT::Grid::extract_grid_t<GridLayerType> GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType<GridType>::v;
  static const bool simplicial_ =
      (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId == Impl::SimplexTopology<dimDomain>::type::id);
  static const bool cubic_ =
      (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId == Impl::CubeTopology<dimDomain>::type::id);
  typedef typename FeMap<GridType, single_geom_, simplicial_, cubic_>::Type FEMapType;

public:
  typedef PDELab::GridFunctionSpace<GridLayerType, FEMapType> BackendType;
  typedef DunePdelabCgMapperWrapper<BackendType> MapperType;
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef BaseFunctionSet::PiolaTransformedDunePdelabWrapper<BackendType,
                                                             EntityType,
                                                             DomainFieldType,
                                                             dimDomain,
                                                             RangeFieldType,
                                                             rangeDim,
                                                             rangeDimCols>
      BaseFunctionSetType;
  static const XT::Grid::Backends layer_backend = XT::Grid::Backends::view;
  static const bool needs_grid_view = true;
  typedef CommunicationChooser<GridLayerType, false> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;

private:
  friend class DunePdelabRtSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class DunePdelabRtSpaceWrapperTraits


} // namespace internal


/**
 * \attention When using alugrid, this space only works on the leaf view, not on level or dd_subdomain* views!
 */
template <class GridViewImp, class RangeFieldImp, size_t rangeDim>
class DunePdelabRtSpaceWrapper<GridViewImp, 0, RangeFieldImp, rangeDim, 1>
    : public RtSpaceInterface<internal::DunePdelabRtSpaceWrapperTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1>,
                              GridViewImp::dimension,
                              rangeDim,
                              1>
{
  typedef DunePdelabRtSpaceWrapper<GridViewImp, 0, RangeFieldImp, rangeDim, 1> ThisType;
  typedef RtSpaceInterface<internal::DunePdelabRtSpaceWrapperTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1>,
                           GridViewImp::dimension,
                           rangeDim,
                           1>
      BaseType;

public:
  typedef internal::DunePdelabRtSpaceWrapperTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1> Traits;

  using BaseType::dimDomain;
  using BaseType::polOrder;

  using typename BaseType::GridLayerType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;

  using typename BaseType::PatternType;
  using typename BaseType::EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  DunePdelabRtSpaceWrapper(GridLayerType grd_vw)
    : grid_view_(grd_vw)
    , fe_map_(grid_view_)
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(Traits::CommunicationChooserType::create(grid_view_))
    , communicator_prepared_(false)
  {
  }

  /**
   * \brief Copy ctor.
   * \note  Manually implemented bc of the std::mutex and our space creation policy
   *        (see https://github.com/pymor/dune-gdt/issues/28)
   */
  DunePdelabRtSpaceWrapper(const ThisType& other)
    : grid_view_(other.grid_view_)
    , fe_map_(grid_view_)
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(Traits::CommunicationChooserType::create(grid_view_))
    , communicator_prepared_(false)
  {
    // make sure our new communicator is prepared if other's was
    if (other.communicator_prepared_)
      const auto& comm DUNE_UNUSED = this->communicator();
  }

  /**
   * \brief Move ctor.
   * \note  Manually implemented bc of the std::mutex and our space creation policy
   *        (see https://github.com/pymor/dune-gdt/issues/28)
   */
  DunePdelabRtSpaceWrapper(ThisType&& source)
    : grid_view_(source.grid_view_)
    , fe_map_(grid_view_)
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(std::move(source.communicator_))
    , communicator_prepared_(source.communicator_prepared_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridLayerType& grid_layer() const
  {
    return grid_view_;
  }

  GridLayerType& grid_layer()
  {
    return grid_view_;
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
      communicator_prepared_ = Traits::CommunicationChooserType::prepare(*this, *communicator_);
    return *communicator_;
  } // ... communicator(...)

private:
  template <class S, size_t d, int p, bool simplicial>
  struct Helper
  {
    static_assert(AlwaysFalse<S>::value, "Not available for this combination!");
  };

  template <class S>
  struct Helper<S, 2, 0, true>
  {
    static std::vector<size_t> local_DoF_indices(const S& space, const EntityType& entity)
    {
      return space.local_DoF_indices_2dsimplex_order0(entity);
    }
  };

public:
  std::vector<size_t> local_DoF_indices(const EntityType& entity) const
  {
    return Helper < ThisType, dimDomain, polOrder,
           Traits::single_geom_ && Traits::simplicial_ > ::local_DoF_indices(*this, entity);
  }

private:
  GridLayerType grid_view_;
  const FEMapType fe_map_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class DunePdelabRtSpaceWrapper< ..., 0, ..., 1 >


#else // HAVE_DUNE_PDELAB


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabRtSpaceWrapper
{
  static_assert(AlwaysFalse<GridViewImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_DUNE_PDELAB_WRAPPER_HH
