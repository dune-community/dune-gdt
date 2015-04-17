// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CG_PDELAB_HH
#define DUNE_GDT_SPACES_CG_PDELAB_HH

#include <memory>

#include <dune/common/fvector.hh>
#include <dune/common/deprecated.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#if HAVE_DUNE_ISTL
# include <dune/istl/paamg/pinfo.hh>
# include <dune/stuff/la/solver/istl_amg.hh>
#endif

#if HAVE_DUNE_PDELAB
# include <dune/pdelab/finiteelementmap/pkfem.hh>
# include <dune/pdelab/finiteelementmap/qkfem.hh>
# include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
# include <dune/pdelab/constraints/conforming.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/spaces/parallel.hh>

#include "../../mapper/pdelab.hh"
#include "../../basefunctionset/pdelab.hh"

#include "interface.hh"

namespace Dune {

namespace GDT {
namespace Spaces {
namespace CG {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class PdelabBased
{
  static_assert(Dune::AlwaysFalse< GridViewImp >::value, "Untested for this combination of dimensions!");
};


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class PdelabBasedTraits
{
public:
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridViewImp                   GridViewType;
  static const int                      polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");
private:
  typedef typename GridViewType::ctype  DomainFieldType;
  static const size_t                   dimDomain = GridViewType::dimension;
public:
  typedef RangeFieldImp                 RangeFieldType;
private:
  template< class G, bool single_geom, bool is_simplex, bool is_cube >
  struct FeMap
  {
    static_assert(Dune::AlwaysFalse< G >::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template< class G >
  struct FeMap< G, true, true, false >
  {
    typedef PDELab::PkLocalFiniteElementMap < GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  template< class G >
  struct FeMap< G, true, false, true >
  {
    typedef PDELab::QkLocalFiniteElementMap < GridViewType, DomainFieldType, RangeFieldType, polOrder> Type;
  };
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                                   == GenericGeometry::SimplexTopology< dimDomain >::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                              == GenericGeometry::CubeTopology< dimDomain >::type::id);
  typedef typename FeMap< GridType, single_geom_, simplicial_, cubic_ >::Type FEMapType;
public:
  typedef PDELab::GridFunctionSpace< GridViewType, FEMapType, PDELab::OverlappingConformingDirichletConstraints > BackendType;
  typedef Mapper::ContinuousPdelabWrapper< BackendType > MapperType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef BaseFunctionSet::PdelabWrapper
      < BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, rangeDim, rangeDimCols >
    BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view = true;

  typedef typename CommunicationChooser<GridViewType>::Type CommunicatorType;
private:
  friend class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >;
}; // class PdelabBasedTraits


template< class GridViewImp, int polynomialOrder, class RangeFieldImp >
class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >
  : public Spaces::CGInterface< PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >,
                                GridViewImp::dimension, 1, 1 >
{
  typedef Spaces::CGInterface< PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >,
                               GridViewImp::dimension, 1, 1 >                    BaseType;
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >       ThisType;
public:
  typedef PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 > Traits;

  static const int    polOrder = Traits::polOrder;
  static const size_t dimDomain = BaseType::dimDomain;
  static const size_t dimRange = BaseType::dimRange;
  static const size_t dimRangeCols = BaseType::dimRangeCols;

  typedef typename Traits::GridViewType         GridViewType;
  typedef typename Traits::RangeFieldType       RangeFieldType;
  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;

  typedef typename GridViewType::ctype              DomainFieldType;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef CommunicationChooser< GridViewType >      CommunicationChooserType;
  typedef typename CommunicationChooserType::Type   CommunicatorType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  typedef typename BaseType::IntersectionType  IntersectionType;
  typedef typename BaseType::EntityType        EntityType;
  typedef typename BaseType::PatternType       PatternType;
  typedef typename BaseType::BoundaryInfoType  BoundaryInfoType;

  explicit PdelabBased(GridViewType gV)
    : gridView_(gV)
    , fe_map_(gridView_)
    , backend_(gridView_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser< GridViewImp >::create(gridView_))
    , communicator_prepared_(false)
  {}

  /**
   * \brief Copy ctor.
   * \note  Manually implemented bc of the std::mutex + communicator_ unique_ptr
   */
  PdelabBased(const ThisType& other)
    : gridView_(other.gridView_)
    , fe_map_(gridView_)
    , backend_(gridView_, fe_map_)
    , mapper_(backend_)
    , communicator_(CommunicationChooser< GridViewImp >::create(gridView_))
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
  PdelabBased(ThisType&& source)
    : gridView_(source.gridView_)
    , fe_map_(source.fe_map_)
    , backend_(source.backend_)
    , mapper_(source.mapper_)
    , communicator_(std::move(source.communicator_))
    , communicator_prepared_(source.communicator_prepared_)
  {}

  ThisType& operator=(const ThisType& other) = delete;

  ThisType& operator=(ThisType&& source) = delete;

  const GridViewType& grid_view() const
  {
    return gridView_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  std::vector< DomainType > lagrange_points(const EntityType& entity) const
  {
    return BaseType::lagrange_points_order_1(entity);
  }

  std::set< size_t > local_dirichlet_DoFs(const EntityType& entity,
                                          const BoundaryInfoType& boundaryInfo) const
  {
    return BaseType::local_dirichlet_DoFs_order_1(entity, boundaryInfo);
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(backend_, entity);
  }

  CommunicatorType& communicator() const
  {
    std::lock_guard< std::mutex > DUNE_UNUSED(gg)(communicator_mutex_);
    if (!communicator_prepared_)
      communicator_prepared_ = CommunicationChooserType::prepare(*this, *communicator_);
    return *communicator_;
  } // ... communicator(...)

private:
  GridViewType gridView_;
  const FEMapType fe_map_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable std::unique_ptr<CommunicatorType> communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class PdelabBased< ..., 1 >


#else // HAVE_DUNE_PDELAB


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class PdelabBased
{
  static_assert(Dune::AlwaysFalse< GridViewImp >::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB


} // namespace CG
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_PDELAB_HH
