// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/constraints/conforming.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/la/container/istl.hh>

#include <dune/gdt/spaces/parallel.hh>

#include <dune/gdt/spaces/basefunctionset/dune-pdelab-wrapper.hh>
#include <dune/gdt/spaces/dg/interface.hh>
#include "../mapper/dune-pdelab-wrapper.hh"

namespace Dune {
namespace GDT {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridViewImp>::value, "Untested for this combination of dimensions!");
};


namespace internal {


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapperTraits
{
public:
  typedef DunePdelabDgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridViewImp GridViewType;
  static const int polOrder    = polynomialOrder;
  static const bool continuous = false;

private:
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

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
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType<GridType>::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                                   == GenericGeometry::SimplexTopology<dimDomain>::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                              == GenericGeometry::CubeTopology<dimDomain>::type::id);
  typedef typename FeMap<GridType, single_geom_, simplicial_, cubic_>::Type FEMapType;

public:
  typedef PDELab::GridFunctionSpace<GridViewType, FEMapType, PDELab::OverlappingConformingDirichletConstraints>
      BackendType;
  typedef Mapper::DiscontinuousPdelabWrapper<BackendType> MapperType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef BaseFunctionSet::DunePdelabWrapper<BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                             rangeDim, rangeDimCols> BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view                       = true;
  typedef CommunicationChooser<GridViewType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;

private:
  friend class DunePdelabDgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class DunePdelabDgSpaceWrapperTraits


} // namespace internal


template <class GridViewImp, int polynomialOrder, class RangeFieldImp>
class DunePdelabDgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>
    : public DGInterface<internal::DunePdelabDgSpaceWrapperTraits<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>,
                         GridViewImp::dimension, 1, 1>
{
  typedef DunePdelabDgSpaceWrapper<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1> ThisType;
  typedef DGInterface<internal::DunePdelabDgSpaceWrapperTraits<GridViewImp, polynomialOrder, RangeFieldImp, 1, 1>,
                      GridViewImp::dimension, 1, 1> BaseType;

public:
  using typename BaseType::Traits;

  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;

private:
  typedef typename Traits::CommunicationChooserType CommunicationChooserType;
  typedef typename Traits::FEMapType FEMapType;

public:
  using typename BaseType::CommunicatorType;

  DunePdelabDgSpaceWrapper(GridViewType gV)
    : grid_view_(gV)
    , fe_map_()
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridViewImp>::create(grid_view_))
    , communicator_prepared_(false)
  {
  }

  /**
   * \brief Copy ctor.
   * \note  Manually implemented bc of the std::mutex + communicator_  unique_ptr
   */
  DunePdelabDgSpaceWrapper(const ThisType& other)
    : grid_view_(other.grid_view_)
    , fe_map_()
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridViewImp>::create(grid_view_))
    , communicator_prepared_(false)
  {
    // make sure our new communicator is prepared if other's was
    if (other.communicator_prepared_)
      const auto& DUNE_UNUSED(comm) = this->communicator();
  }

  /**
   * \brief Move ctor.
   * \note  Manually implemented bc of the std::mutex.
   */
  DunePdelabDgSpaceWrapper(ThisType&& source)
    : grid_view_(source.grid_view_)
    , fe_map_(source.fe_map_)
    , backend_(source.backend_)
    , mapper_(source.mapper_)
    , communicator_(std::move(source.communicator_))
    , communicator_prepared_(source.communicator_prepared_)
  {
  }

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridViewType& grid_view() const
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
    std::lock_guard<std::mutex> DUNE_UNUSED(gg)(communicator_mutex_);
    if (!communicator_prepared_)
      communicator_prepared_ = CommunicationChooser<GridViewType>::prepare(*this, *communicator_);
    return *communicator_;
  } // ... communicator(...)

private:
  GridViewType grid_view_;
  const FEMapType fe_map_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class DunePdelabDgSpaceWrapper< ..., 1 >


#else // HAVE_DUNE_PDELAB


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class DunePdelabDgSpaceWrapper
{
  static_assert(Dune::AlwaysFalse<GridViewImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_DG_DUNE_PDELAB_WRAPPER_HH
