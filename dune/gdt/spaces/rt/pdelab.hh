// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_RT_PDELAB_HH
#define DUNE_GDT_SPACES_RT_PDELAB_HH

#include <type_traits>
#include <limits>
#include <mutex>

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/basefunctionset/pdelab.hh>
#include <dune/gdt/mapper/pdelab.hh>
#include <dune/gdt/spaces/parallel.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace RT {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class PdelabBased
{
  static_assert(AlwaysFalse<GridViewImp>::value, "Untested for these dimensions or polynomial order!");
}; // class PdelabBased


namespace internal {


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class PdelabBasedTraits
{
public:
  typedef PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols> derived_type;
  typedef GridViewImp GridViewType;
  static const int polOrder    = polynomialOrder;
  static const bool continuous = false;
  static_assert(polOrder == 0, "Untested!");
  static_assert(rangeDim == GridViewType::dimension, "Untested!");
  static_assert(rangeDimCols == 1, "Untested!");

private:
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

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
    typedef PDELab::RaviartThomasLocalFiniteElementMap<GridViewType, DomainFieldType, RangeFieldType, polOrder,
                                                       Dune::GeometryType::simplex> Type;
  };
  template <class G>
  struct FeMap<G, true, false, true>
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap<GridViewType, DomainFieldType, RangeFieldType, polOrder,
                                                       Dune::GeometryType::cube> Type;
  };
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType<GridType>::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                                   == GenericGeometry::SimplexTopology<dimDomain>::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType<GridType>::topologyId
                              == GenericGeometry::CubeTopology<dimDomain>::type::id);
  typedef typename FeMap<GridType, single_geom_, simplicial_, cubic_>::Type FEMapType;

public:
  typedef PDELab::GridFunctionSpace<GridViewType, FEMapType> BackendType;
  typedef Mapper::ContinuousPdelabWrapper<BackendType> MapperType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef BaseFunctionSet::PiolaTransformedPdelabWrapper<BackendType, EntityType, DomainFieldType, dimDomain,
                                                         RangeFieldType, rangeDim, rangeDimCols> BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view                       = true;
  typedef CommunicationChooser<GridViewType> CommunicationChooserType;
  typedef typename CommunicationChooserType::Type CommunicatorType;

private:
  friend class PdelabBased<GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols>;
}; // class PdelabBasedTraits


} // namespace internal


template <class GridViewImp, class RangeFieldImp, size_t rangeDim>
class PdelabBased<GridViewImp, 0, RangeFieldImp, rangeDim, 1>
    : public Spaces::RTInterface<internal::PdelabBasedTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1>,
                                 GridViewImp::dimension, rangeDim, 1>
{
  typedef PdelabBased<GridViewImp, 0, RangeFieldImp, rangeDim, 1> ThisType;
  typedef Spaces::RTInterface<internal::PdelabBasedTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1>,
                              GridViewImp::dimension, rangeDim, 1> BaseType;

public:
  typedef internal::PdelabBasedTraits<GridViewImp, 0, RangeFieldImp, rangeDim, 1> Traits;

  using BaseType::dimDomain;
  using BaseType::polOrder;

  using typename BaseType::GridViewType;
  using typename BaseType::BackendType;
  using typename BaseType::MapperType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::CommunicatorType;

  using typename BaseType::PatternType;
  using typename BaseType::EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  PdelabBased(GridViewType gV)
    : grid_view_(gV)
    , fe_map_(grid_view_)
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridViewType>::create(grid_view_))
    , communicator_prepared_(false)
  {
  }

  /**
   * \brief Copy ctor.
   * \note  Manually implemented bc of the std::mutex and our space creation policy
   *        (see https://github.com/pymor/dune-gdt/issues/28)
   */
  PdelabBased(const ThisType& other)
    : grid_view_(other.grid_view_)
    , fe_map_(grid_view_)
    , backend_(grid_view_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser<GridViewType>::create(grid_view_))
    , communicator_prepared_(false)
  {
    // make sure our new communicator is prepared if other's was
    if (other.communicator_prepared_)
      const auto& DUNE_UNUSED(comm) = this->communicator();
  }

  /**
   * \brief Move ctor.
   * \note  Manually implemented bc of the std::mutex and our space creation policy
   *        (see https://github.com/pymor/dune-gdt/issues/28)
   */
  PdelabBased(ThisType&& source)
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
  GridViewType grid_view_;
  const FEMapType fe_map_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class PdelabBased< ..., 0, ..., 1 >


#else // HAVE_DUNE_PDELAB


template <class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1>
class PdelabBased
{
  static_assert(AlwaysFalse<GridViewImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB


} // namespace RT
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_PDELAB_HH
