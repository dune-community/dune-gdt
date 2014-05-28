// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH
#define DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH

#include <type_traits>
#include <limits>

#include <dune/common/static_assert.hh>

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/capabilities.hh>

#if HAVE_DUNE_PDELAB
# include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
# include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#endif // HAVE_DUNE_PDELAB

#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/exceptions.hh>

#include <dune/gdt/basefunctionset/pdelab.hh>
#include <dune/gdt/mapper/pdelab.hh>

#include "../../../spaces/interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {
namespace RaviartThomas {

#if HAVE_DUNE_PDELAB


// forward, to be used in the traits and to allow for specialization
template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert(AlwaysFalse< GridViewImp >::value, "Untested for these dimensions or polynomial order!");
}; // class PdelabBased


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols >
class PdelabBasedTraits
{
public:
  typedef PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef GridViewImp GridViewType;
  static const int    polOrder = polynomialOrder;
  static_assert(polOrder == 0, "Untested!");
  static_assert(rangeDim == GridViewType::dimension, "Untested!");
  static_assert(rangeDimCols == 1, "Untested!");
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
    static_assert(AlwaysFalse< G >::value,
                  "This space is only implemented for either fully simplicial or fully cubic grids!");
  };
  template< class G >
  struct FeMap< G, true, true, false >
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap< GridViewType, DomainFieldType,
                                                        RangeFieldType, polOrder, Dune::GeometryType::simplex > Type;
  };
  template< class G >
  struct FeMap< G, true, false, true >
  {
    typedef PDELab::RaviartThomasLocalFiniteElementMap< GridViewType, DomainFieldType,
                                                        RangeFieldType, polOrder, Dune::GeometryType::cube > Type;
  };
  typedef typename GridViewType::Grid GridType;
  static const bool single_geom_ = Dune::Capabilities::hasSingleGeometryType< GridType >::v;
  static const bool simplicial_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                                   == GenericGeometry::SimplexTopology< dimDomain >::type::id);
  static const bool cubic_ = (Dune::Capabilities::hasSingleGeometryType< GridType >::topologyId
                              == GenericGeometry::CubeTopology< dimDomain >::type::id);
  typedef typename FeMap< GridType, single_geom_, simplicial_, cubic_ >::Type FEMapType;
public:
  typedef PDELab::GridFunctionSpace< GridViewType, FEMapType > BackendType;
  typedef Mapper::SimplePdelabWrapper< BackendType > MapperType;
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef BaseFunctionSet::PiolaTransformedPdelabWrapper
      < BackendType, EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >
    BaseFunctionSetType;
  static const Stuff::Grid::ChoosePartView part_view_type = Stuff::Grid::ChoosePartView::view;
  static const bool needs_grid_view = true;
  typedef double CommunicatorType;
private:
  friend class PdelabBased< GridViewImp, polynomialOrder, RangeFieldImp, rangeDim, rangeDimCols >;
}; // class PdelabBasedTraits


template< class GridViewImp, class RangeFieldImp, int rangeDim >
class PdelabBased< GridViewImp, 0, RangeFieldImp, rangeDim, 1 >
  : public SpaceInterface< PdelabBasedTraits< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > >
{
  typedef PdelabBased< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > ThisType;
public:
  typedef PdelabBasedTraits< GridViewImp, 0, RangeFieldImp, rangeDim, 1 > Traits;

  typedef typename Traits::GridViewType GridViewType;
  static const int                      polOrder = Traits::polOrder;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridViewType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;

private:
  typedef typename Traits::FEMapType FEMapType;

public:
  PdelabBased(const std::shared_ptr< const GridViewType >& gV)
    : grid_view_(gV)
    , fe_map_(std::make_shared< FEMapType >(*(grid_view_)))
    , backend_(std::make_shared< BackendType >(const_cast< GridViewType& >(*grid_view_), *fe_map_))
    , mapper_(std::make_shared< MapperType >(*backend_))
    , communicator_(0.0)
  {}

  PdelabBased(const ThisType& other)
    : grid_view_(other.grid_view_)
    , fe_map_(other.fe_map_)
    , backend_(other.backend_)
    , mapper_(other.mapper_)
    , communicator_(0.0)
  {}

  ThisType& operator=(const ThisType& other)
  {
    if (this != &other) {
      grid_view_ = other.grid_view_;
      fe_map_ = other.fe_map_;
      backend_ = other.backend_;
      mapper_ = other.mapper_;
    }
    return *this;
  }

  ~PdelabBased() {}

  const std::shared_ptr< const GridViewType >& grid_view() const
  {
    return grid_view_;
  }

  const BackendType& backend() const
  {
    return *backend_;
  }

  const MapperType& mapper() const
  {
    return *mapper_;
  }

  BaseFunctionSetType base_function_set(const EntityType& entity) const
  {
    return BaseFunctionSetType(*backend_, entity);
  }

  /**
   *  \brief  Computes a vector 'indices' of length entity.count< 1 >(), where 'indices[intersection.indexInInside()]'
   *          is the index of the basis function (aka the local DoF index) corresponding to the intersection.
   */
  std::vector< size_t > local_DoF_indices(const EntityType& entity) const
  {
    static_assert(dimDomain == 2 && polOrder == 0, "Not implemented!");
    // prepare
    const size_t num_intersections = entity.template count< 1 >();
    std::vector< size_t > local_DoF_index_of_vertex(num_intersections, std::numeric_limits< size_t >::infinity());
    std::vector< size_t > local_DoF_index_of_intersection(num_intersections, std::numeric_limits< size_t >::infinity());
    typedef typename BaseFunctionSetType::DomainType DomainType;
    std::vector< DomainType > vertices(num_intersections, DomainType(0));
    std::vector< bool > lies_on_intersection(num_intersections, false);
    DomainType corner(0);
    typedef typename BaseFunctionSetType::RangeType RangeType;
    const RangeType one(1);
    std::vector< RangeType > basis_values(num_intersections, one);
    const auto basis = base_function_set(entity);
    assert(basis.size() == num_intersections);
    const auto geometry = entity.geometry();
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(geometry.type());
    // find the basis function index that corresponds to each vertex of the entity
    // (find the basis function that evaluates to zero at the vertex, and nonzero at the other ones)
    // therefore we walk the vertices
    assert(int(num_intersections) == entity.template count< dimDomain >());
    for (size_t vv = 0; vv < num_intersections; ++vv) {
      const auto vertex_ptr = entity.template subEntity< dimDomain >(int(vv));
      const auto& vertex = *vertex_ptr;
      // get the vertex coordinates
      vertices[vv] = vertex.geometry().center();
      const auto& vertex_entity = reference_element.position(int(vv), dimDomain);
      // evalaute the basis
      basis.evaluate(vertex_entity, basis_values);
      // and find the basis that evaluates zero here
      size_t zeros = 0;
      size_t nonzeros = 0;
      for (size_t ii = 0; ii < num_intersections; ++ii) {
        // we would like to check against 0, but there is a bug in dune-commons FloatCmp
        if (Stuff::Common::FloatCmp::eq(basis_values[ii] + one, one)) {
          // this is a candidate for the basis function we are looking for
          local_DoF_index_of_vertex[vv] = ii;
          ++zeros;
        } else
          ++nonzeros;
      }
      // make sure there was only one candidate
      if (zeros != 1 || nonzeros != (num_intersections - 1))
        DUNE_THROW_COLORFULLY(Stuff::Exceptions::internal_error,
                              "This must not happen for RTN0 in 2d!\n"
                              << "  zeros    = " << zeros << "\n"
                              << "  nonzeros = " << nonzeros);
    } // walk the vertices
    // so from here on we have the local DoF index that corresponds to each vertex vv in local_DoF_index_of_vertex[vv]
    // now we need to find the intersection opposite to this vertex
    // therefore we walk the intersections
    size_t intersection_counter = 0;
    const auto intersection_it_end = grid_view_->iend(entity);
    for (auto intersection_it = grid_view_->ibegin(entity);
         intersection_it != intersection_it_end;
         ++intersection_it) {
      const auto& intersection = *intersection_it;
      const auto intersection_geometry = intersection.geometry();
      const size_t local_intersection_index = intersection.indexInInside();
      // make sure this index has not been already taken by another intersection
      assert(local_DoF_index_of_intersection[local_intersection_index] == std::numeric_limits< size_t >::infinity());
      // walk the corners of the intersection
      for (size_t cc = 0; cc < num_intersections; ++cc) {
        corner = intersection_geometry.corner(int(cc));
        // check which vertices lie on the intersection
        for (size_t vv = 0; vv < num_intersections; ++vv)
          if (Stuff::Common::FloatCmp::eq(vertices[vv], corner))
            lies_on_intersection[vv] = true;
      } // walk the corners of the intersection
      // now see if we find a vertex that does not lie on the intersection
      size_t found = 0;
      size_t missed = 0;
      for (size_t vv = 0; vv < num_intersections; ++vv) {
        if (!(lies_on_intersection[vv])) {
          // this is a good candidate
          // so the local DoF id that corresponds to this vertex is the one that corresponds to the intersection
          local_DoF_index_of_intersection[local_intersection_index] = local_DoF_index_of_vertex[vv];
          ++found;
        } else
          ++missed;
        // and clear for the next intersection
        lies_on_intersection[vv] = false;
      } // walk the vertices of this entity
      // make sure there was only one candidate
      if (found != 1 || missed != (num_intersections - 1))
        DUNE_THROW_COLORFULLY(Stuff::Exceptions::internal_error,
                              "This must not happen for RTN0 in 2d!\n"
                              << "  found  = " << found << "\n"
                              << "  missed = " << missed);
      ++intersection_counter;
    } // walk the intersection
    assert(intersection_counter == num_intersections);
    return local_DoF_index_of_intersection;
  } // ... local_DoF_indices(...)

  double& communicator() const
  {
    return communicator_;
  }

private:
  std::shared_ptr< const GridViewType > grid_view_;
  std::shared_ptr< const FEMapType > fe_map_;
  std::shared_ptr< const BackendType > backend_;
  std::shared_ptr< const MapperType > mapper_;
  mutable double communicator_;
}; // class PdelabBased< ..., 0, ..., 1 >


#else // HAVE_DUNE_PDELAB


template< class GridViewImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class PdelabBased
{
  static_assert(AlwaysFalse< GridViewImp >::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace RaviartThomas
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_RAVIARTTHOMASSPACE_PDELAB_HH
