// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_CONTINUOUS_LAGRANGE_FEM_LOCALFUNCTIONS_HH
#define DUNE_GDT_SPACE_CONTINUOUS_LAGRANGE_FEM_LOCALFUNCTIONS_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>

#include <dune/localfunctions/lagrange/equidistantpoints.hh>
#include <dune/localfunctions/lagrange.hh>

#include <dune/fem/space/common/allgeomtypes.hh>

#include <dune/fem_localfunctions/localfunctions/transformations.hh>
#include <dune/fem_localfunctions/basefunctions/genericbasefunctionsetstorage.hh>
#include <dune/fem_localfunctions/basefunctionsetmap/basefunctionsetmap.hh>
#include <dune/fem_localfunctions/space/genericdiscretefunctionspace.hh>

#include <dune/stuff/common/color.hh>
#include <dune/stuff/common/float_cmp.hh>

#include "../../mapper/fem.hh"
#include "../../basefunctionset/fem-localfunctions.hh"
#include "../constraints.hh"
#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace ContinuousLagrangeSpace {


// forward, to be used in the traits and to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapper;


// forward, to allow for specialization
template< class GridPartImp, int polynomialOrder, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemLocalfunctionsWrapperTraits;


/**
 *  \brief Traits class for ContinuousLagrangeSpace for dimRange 1x1.
 */
template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
{
public:
  typedef GridPartImp                   GridPartType;
  static const int                      polOrder = polynomialOrder;
  dune_static_assert((polOrder >= 1), "ERROR: wrong polOrder given!");
private:
  typedef typename GridPartType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridPartType::dimension;
public:
  typedef RangeFieldImp                 RangeFieldType;
  static const unsigned int             dimRange = 1;
  static const unsigned int             dimRangeCols = 1;
  typedef FemLocalfunctionsWrapper< GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols > derived_type;
private:
  typedef FemLocalfunctionsWrapperTraits< GridPartType, polOrder, RangeFieldType, dimRange, dimRangeCols > ThisType;
public:
  typedef Dune::LagrangeLocalFiniteElement< Dune::EquidistantPointSet,
                                            dimDomain,
                                            DomainFieldType,
                                            RangeFieldType >        FiniteElementType;
private:
  typedef Dune::FemLocalFunctions::BaseFunctionSetMap<  GridPartType,
                                                        FiniteElementType,
                                                        Dune::FemLocalFunctions::NoTransformation,
                                                        Dune::FemLocalFunctions::SimpleStorage,
                                                        polOrder,
                                                        polOrder >  BaseFunctionSetMapType;
public:
  typedef Dune::FemLocalFunctions::DiscreteFunctionSpace< BaseFunctionSetMapType >  BackendType;
  typedef Mapper::FemDofWrapper< typename BackendType::MapperType >                 MapperType;
  typedef BaseFunctionSet::FemLocalfunctionsWrapper< BaseFunctionSetMapType,
              DomainFieldType, dimDomain, RangeFieldType, dimRange, dimRangeCols >  BaseFunctionSetType;
  typedef typename BaseFunctionSetType::EntityType                                  EntityType;
private:
  template< class G, int p, class R, int r, int rC >
  friend class FemLocalfunctionsWrapper;
}; // class FemLocalfunctionsWrapperTraits< ..., 1, 1 >


template< class GridPartImp, int polynomialOrder, class RangeFieldImp >
class FemLocalfunctionsWrapper< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 >
  : public SpaceInterface< FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > >
{
  typedef SpaceInterface< FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > > BaseType;
public:
  typedef FemLocalfunctionsWrapperTraits< GridPartImp, polynomialOrder, RangeFieldImp, 1, 1 > Traits;

  typedef typename Traits::GridPartType   GridPartType;
  typedef typename GridPartType::ctype    DomainFieldType;
  static const int                        polOrder = Traits::polOrder;
  static const unsigned int               dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const unsigned int               dimRangeCols = Traits::dimRangeCols;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  typedef Dune::Stuff::LA::SparsityPatternDefault PatternType;

private:
  typedef typename Traits::BaseFunctionSetMapType BaseFunctionSetMapType;

public:
  FemLocalfunctionsWrapper(const GridPartType& gridP)
    : gridPart_(assertGridPart(gridP))
    , baseFunctionSetMap_(gridPart_)
    , backend_(const_cast< GridPartType& >(gridPart_), baseFunctionSetMap_)
    , mapper_(backend_.mapper())
    , tmp_global_indices_(mapper_.maxNumDofs())
  {}

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  bool continuous() const
  {
    return true;
  }

  const MapperType& mapper() const
  {
    return mapper_;
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return BaseFunctionSetType(baseFunctionSetMap_, entity);
  }

  template< class R >
  void localConstraints(const EntityType& /*entity*/,
                        Constraints::LocalDefault< R >& /*ret*/) const
  {
    dune_static_assert(Dune::AlwaysFalse< R >::value, "ERROR: not implemented for arbitrary constraints!");
  }

  void localConstraints(const EntityType& entity,
                        Constraints::Dirichlet< typename GridPartType::IntersectionType, RangeFieldType, true >& ret) const
  {
    static_assert(dimDomain == 2, "This does not work for other dimensions!");
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    // only work if we are at the boundary
    if (entity.hasBoundaryIntersections()) {
      // get local finite elements
      // * we are a CG space
      const auto cg_finite_element = backend_.finiteElement(entity);
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      // * but we also need a local DG finite element
      typedef Dune::DGLocalFiniteElement< typename Traits::FiniteElementType > DgFiniteElementType;
      const DgFiniteElementType dg_finite_element(entity.type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");

      // first we loop over all vertices of the entity
      typedef FieldVector< DomainFieldType, dimDomain > DomainType;
      std::vector< DomainType > global_vertices(entity.template count< dimDomain >(), DomainType(0));
      std::vector< size_t > local_DoF_ids(global_vertices.size(), 0);
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get the vertex
        const auto vertexPtr = entity.template subEntity< dimDomain >(local_vertex_id);
        const auto& vertex = *vertexPtr;
        global_vertices[local_vertex_id] = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            // get the local DoF to this vertex
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            local_DoF_ids[local_vertex_id] = entity_dg_local_key.index();
            // there must be one and only one for a polorder 1 lagrange basis
            break;
          }
        } // loop over all local DoFs
      } // loop over all vertices of the entity

      // then we walk the intersections
      std::set< size_t > local_dirichlet_DoF_ids;
      const auto intersection_it_end = gridPart_.iend(entity);
      for (auto intersection_it = gridPart_.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (ret.gridBoundary().dirichlet(intersection)) {
          const auto& intersection_geometry = intersection.geometry();
          // and walk its corners (i.e. the vertices in 2d)
          for (size_t local_intersection_corner_id = 0;
               int(local_intersection_corner_id) < intersection_geometry.corners();
               ++local_intersection_corner_id) {
            const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
            // to check which vertex this corner is
            // loop over all vertices of the entity again
            for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
              // and check for equality
              if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id])) {
                // this vertex is on the dirichlet boundary, so we add the local DoF id this vertex corresponds to
                local_dirichlet_DoF_ids.insert(local_DoF_ids[local_vertex_id]);
              }
            } // loop over all vertices of the entity
          } // walk its corners
        }
      } // walk the intersections

      // finally we compute the local constraints
      const size_t num_rows = local_dirichlet_DoF_ids.size();
      if (num_rows > 0) {
        const size_t num_cols = baseFunctionSet(entity).size();
        mapper_.globalIndices(entity, tmp_global_indices_);
        ret.setSize(num_rows, num_cols);
        size_t local_row = 0;
        const RangeFieldType zero(0);
        const RangeFieldType one(1);
        for (const size_t& local_dirichlet_DoF_id : local_dirichlet_DoF_ids) {
          ret.globalRow(local_row) = tmp_global_indices_[local_dirichlet_DoF_id];
          for (size_t jj = 0; jj < ret.cols(); ++jj) {
            ret.globalCol(jj) = tmp_global_indices_[jj];
            if (tmp_global_indices_[jj] == tmp_global_indices_[local_dirichlet_DoF_id])
              ret.value(local_row, jj) = one;
            else
              ret.value(local_row, jj) = zero;
          }
          ++local_row;
        }
      } else {
        ret.setSize(0, 0);
      }
    }
  } // ... localConstraints(...)

  void localConstraints(const EntityType& entity,
                        Constraints::Dirichlet< typename GridPartType::IntersectionType, RangeFieldType, false >& ret) const
  {
    static_assert(dimDomain == 2, "This does not work for other dimensions!");
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    // only work if we are at the boundary
    if (entity.hasBoundaryIntersections()) {
      // get local finite elements
      // * we are a CG space
      const auto cg_finite_element = backend_.finiteElement(entity);
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      // * but we also need a local DG finite element
      typedef Dune::DGLocalFiniteElement< typename Traits::FiniteElementType > DgFiniteElementType;
      const DgFiniteElementType dg_finite_element(entity.type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");

      // first we loop over all vertices of the entity
      typedef FieldVector< DomainFieldType, dimDomain > DomainType;
      std::vector< DomainType > global_vertices(entity.template count< dimDomain >(), DomainType(0));
      std::vector< size_t > local_DoF_ids(global_vertices.size(), 0);
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get the vertex
        const auto vertexPtr = entity.template subEntity< dimDomain >(local_vertex_id);
        const auto& vertex = *vertexPtr;
        global_vertices[local_vertex_id] = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(ii);
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            // get the local DoF to this vertex
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(ii);
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            local_DoF_ids[local_vertex_id] = entity_dg_local_key.index();
            // there must be one and only one for a polorder 1 lagrange basis
            break;
          }
        } // loop over all local DoFs
      } // loop over all vertices of the entity

      // then we walk the intersections
      std::set< size_t > local_dirichlet_DoF_ids;
      const auto intersection_it_end = gridPart_.iend(entity);
      for (auto intersection_it = gridPart_.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (ret.gridBoundary().dirichlet(intersection)) {
          const auto& intersection_geometry = intersection.geometry();
          // and walk its corners (i.e. the vertices in 2d)
          for (size_t local_intersection_corner_id = 0;
               int(local_intersection_corner_id) < intersection_geometry.corners();
               ++local_intersection_corner_id) {
            const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
            // to check which vertex this corner is
            // loop over all vertices of the entity again
            for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
              // and check for equality
              if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id])) {
                // this vertex is on the dirichlet boundary, so we add the local DoF id this vertex corresponds to
                local_dirichlet_DoF_ids.insert(local_DoF_ids[local_vertex_id]);
              }
            } // loop over all vertices of the entity
          } // walk its corners
        }
      } // walk the intersections

      // finally we compute the local constraints
      const size_t num_rows = local_dirichlet_DoF_ids.size();
      if (num_rows > 0) {
        const size_t num_cols = baseFunctionSet(entity).size();
        mapper_.globalIndices(entity, tmp_global_indices_);
        ret.setSize(num_rows, num_cols);
        size_t local_row = 0;
        const RangeFieldType zero(0);
        for (const size_t& local_dirichlet_DoF_id : local_dirichlet_DoF_ids) {
          ret.globalRow(local_row) = tmp_global_indices_[local_dirichlet_DoF_id];
          for (size_t jj = 0; jj < ret.cols(); ++jj) {
            ret.globalCol(jj) = tmp_global_indices_[jj];
            ret.value(local_row, jj) = zero;
          }
          ++local_row;
        }
      } else {
        ret.setSize(0, 0);
      }
    }
  } // ... localConstraints(...)

  using BaseType::computePattern;

  template< class LocalGridPartType, class OtherSpaceType >
  PatternType* computePattern(const LocalGridPartType& localGridPart,
                              const OtherSpaceType& otherSpace) const
  {
    return BaseType::computeVolumePattern(localGridPart, otherSpace);
  }

private:
  static const GridPartType& assertGridPart(const GridPartType& gP)
  {
    // check
    typedef typename Dune::Fem::AllGeomTypes< typename GridPartType::IndexSetType,
                                              typename GridPartType::GridType > AllGeometryTypes;
    const AllGeometryTypes allGeometryTypes(gP.indexSet());
    const std::vector< Dune::GeometryType >& geometryTypes = allGeometryTypes.geomTypes(0);
    if (!(geometryTypes.size() == 1 && geometryTypes[0].isSimplex()))
      DUNE_THROW(Dune::NotImplemented,
                 "\n" << Dune::Stuff::Common::colorStringRed("ERROR:")
                 << " this space is only implemented for simplicial grids!");
    return gP;
  } // ... assertGridPart(...)

  const GridPartType& gridPart_;
  BaseFunctionSetMapType baseFunctionSetMap_;
  const BackendType backend_;
  const MapperType mapper_;
  mutable Dune::DynamicVector< size_t > tmp_global_indices_;
}; // class FemLocalfunctionsWrapper< ..., 1, 1 >


} // namespace ContinuousLagrangeSpace
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_CONTINUOUS_LAGRANGE_FEM_LOCALFUNCTIONS_HH
