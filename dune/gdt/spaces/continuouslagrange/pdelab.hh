// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH

#include <memory>
#include <type_traits>

#include <dune/geometry/genericgeometry/topologytypes.hh>

#include <dune/common/typetraits.hh>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_ISTL
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/istl/paamg/pinfo.hh>
# include <dune/stuff/common/reenable_warnings.hh>
# include <dune/stuff/la/solver/istl_amg.hh>
#endif

#if HAVE_DUNE_PDELAB
# include <dune/stuff/common/disable_warnings.hh>
#   include <dune/pdelab/finiteelementmap/pkfem.hh>
#   include <dune/pdelab/finiteelementmap/qkfem.hh>
#   include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#   include <dune/pdelab/constraints/conforming.hh>
# include <dune/stuff/common/reenable_warnings.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/gdt/spaces/parallel.hh>

#include "../../mapper/pdelab.hh"
#include "../../basefunctionset/pdelab.hh"

#include "base.hh"

namespace Dune {

namespace GDT {
namespace Spaces {
namespace ContinuousLagrange {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert((Dune::AlwaysFalse< GridViewImp >::value), "Untested for this combination of dimensions!");
};


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBasedTraits
{
public:
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridViewImp                   GridViewType;
  static const int                      polOrder = polynomialOrder;
  static_assert(polOrder >= 1, "Wrong polOrder given!");
private:
  typedef typename GridViewType::ctype  DomainFieldType;
public:
  static const unsigned int             dimDomain = GridViewType::dimension;
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = rangeDim;
  static const unsigned int             dimRangeCols = rangeDimCols;
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
  typedef Mapper::SimplePdelabWrapper< BackendType > MapperType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef BaseFunctionSet::PdelabWrapper
      < BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >
    BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view = true;

  typedef typename CommunicationChooser<GridViewType>::Type CommunicatorType;
private:
  friend class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >;
}; // class SpaceWrappedFemContinuousLagrangeTraits


template< class GridViewImp, int polynomialOrder, class RangeFieldImp >
class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >
  : public Spaces::ContinuousLagrangeBase< PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >
                                      , GridViewImp::dimension, RangeFieldImp, 1, 1 >
{
  typedef Spaces::ContinuousLagrangeBase< PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >
                                     , GridViewImp::dimension, RangeFieldImp, 1, 1 >  BaseType;
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 >          ThisType;
public:
  typedef PdelabBasedTraits< GridViewImp, polynomialOrder, RangeFieldImp, 1, 1 > Traits;

  typedef typename Traits::GridViewType GridViewType;
  static const int                      polOrder = Traits::polOrder;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridViewType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::CommunicatorType     CommunicatorType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  typedef typename BaseType::IntersectionType  IntersectionType;
  typedef typename BaseType::EntityType        EntityType;
  typedef typename BaseType::PatternType       PatternType;
  typedef typename BaseType::BoundaryInfoType  BoundaryInfoType;

  PdelabBased(const std::shared_ptr< const GridViewType >& gV)
    : gridView_(gV)
    , fe_map_(std::make_shared< FEMapType >(*(gridView_)))
    , backend_(std::make_shared< BackendType >(const_cast< GridViewType& >(*gridView_), *(*fe_map_)))
    , mapper_(std::make_shared< MapperType >(*(*backend_)))
    , communicator_(CommunicationChooser<GridViewImp>::create(*gridView_))
    , communicator_prepared_(false)
  {}

  PdelabBased(const ThisType& other)
    : gridView_(other.gridView_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , communicator_(other.communicator_)
    , communicator_prepared_(other.communicator_prepared_)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      gridView_ = other.gridView_;
      fe_map_ = other.fe_map_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
      communicator_ = other.communicator_;
      communicator_prepared_ = other.communicator_prepared_;
    }
    return *this;
  }

  ~PdelabBased() {}

  const std::shared_ptr< const GridViewType >& grid_view() const
  {
    return gridView_;
  }

  const BackendType& backend() const
  {
    return *(*backend_);
  }

  const MapperType& mapper() const
  {
    return *(*mapper_);
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(*(*backend_), entity);
  }

  CommunicatorType& communicator() const
  {
    std::lock_guard<std::mutex> gg(communicator_mutex_);
    if (!communicator_prepared_) {
      communicator_prepared_ = CommunicationChooser<GridViewType>::prepare(*this, *communicator_);
    }
    return *communicator_;
  } // ... communicator(...)

private:
  std::shared_ptr< const GridViewType > gridView_;
  DS::PerThreadValue<std::shared_ptr< const FEMapType >> fe_map_;
  DS::PerThreadValue<std::shared_ptr< const BackendType >> backend_;
  DS::PerThreadValue<std::shared_ptr< const MapperType >> mapper_;
  mutable std::shared_ptr< CommunicatorType > communicator_;
  mutable bool communicator_prepared_;
  mutable std::mutex communicator_mutex_;
}; // class PdelabBased< ..., 1 >


#else // HAVE_DUNE_PDELAB


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert((Dune::AlwaysFalse< GridViewImp >::value), "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace ContinuousLagrange
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_PDELAB_HH
