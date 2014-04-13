// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_OSWALD_HH
#define DUNE_GDT_OPERATORS_OSWALD_HH

#include <vector>
#include <set>
#include <limits>

#include <dune/stuff/common/vector.hh>

#include <dune/gdt/assembler/gridwalker.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/discontinuouslagrange/fem-localfunctions.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


#if 0
// forward, to be used in the traits
template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp = double >
class OswaldInterpolationLocalizable;


template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class OswaldInterpolationLocalizableTraits
{
public:
  typedef OswaldInterpolationLocalizable< GridViewImp, SourceImp, RangeImp, FieldImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef SourceImp   SourceType;
  typedef RangeImp    RangeType;
  typedef FieldImp    FieldType;
private:
  static_assert(std::is_base_of< ConstDiscreteFunction< typename SourceType::SpaceType
                                                      , typename SourceType::VectorType >
                               , RangeType >::value,
                "SourceType has to be a ConstDiscreteFunction!");
  typedef typename SourceType::SpaceType S;
  static_assert(std::is_same< S
                            , DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< typename S::GridPartType
                                                                                  , 1
                                                                                  , FieldType, 1, 1 > >::value,
                "The SpaceType of SourceType has to be a DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper!");
  typedef typename RangeType::SpaceType R;
  static_assert(std::is_same< R
                            , DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< typename R::GridPartType
                                                                                  , 1
                                                                                  , FieldType, 1, 1 > >::value,
                "The SpaceType of RangeType has to be a DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper!");
}; // class OswaldInterpolationLocalizableTraits


template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class OswaldInterpolationLocalizable
  : public LocalizableOperatorInterface< OswaldInterpolationLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > >
  , public Functor::Codim0And1< GridViewImp >
{
  typedef Functor::Codim0And1< GridViewImp > FunctorType;
public:
  typedef OswaldInterpolationLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > Traits;

  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::FieldType    FieldType;

  using typename FunctorType::EntityType;
  using typename FunctorType::IntersectionType;

private:
  typedef typename RangeType::SpaceType SpaceType;
  static const unsigned int               dimDomain = RangeType::dimDomain;
  typedef typename RangeType::DomainType  DomainType;

  typedef typename SpaceType::Traits::ContinuousFiniteElementType CgFemType;
  typedef typename SpaceType::Traits::FiniteElementType           DgFemType;

public:
  OswaldInterpolationLocalizable(const GridViewType& grid_view, const SourceType& source, RangeType& range)
    : grid_view_(grid_view)
    , source_(source)
    , range_(range)
    , prepared_(false)
    , applied_(false)
    , finalized_(false)
    , current_entity_(std::numeric_limits< size_t >::max())
    , cg_fem_(nullptr)
    , dg_fem_(nullptr)
  {
    assert(source.space().mapper().size() == range.space().mapper().size());
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

  void apply()
  {
    if (!applied_) {
      GridWalker< GridViewType > grid_walker(grid_view_);
      grid_walker.add(*this, new ApplyOn::BoundaryIntersections< GridViewType >());
      grid_walker.walk();
      applied_ = true;
    }
  } // ... apply()

  virtual void prepare() /*DS_OVERRIDE DS_FINAL*/
  {
    if (!prepared_) {
      num_local_DoFs_per_global_DoF_ = std::vector< size_t >(range_.space().mapper().size(), 0);
      vertex_to_dof_id_map_ = std::vector< std::set< size_t > >(grid_view_.indexSet().size(dimDomain));
      Stuff::Common::clear(range_.vector());
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) /*DS_OVERRIDE DS_FINAL*/
  {
    assert(prepared_);
    if (!finalized_) {
      // some preparations
      prepare_for_current_(entity);
      const auto local_source = source_.local_discrete_function(entity);
      const auto local_source_DoF_vector = local_source.vector();
      const size_t num_vertices = entity.template count< dimDomain >();
      assert(cg_fem_);
      assert(dg_fem_);
      const auto& dg_local_coefficients = dg_fem_->localCoefficients();
      const auto& cg_local_coefficients = cg_fem_->localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // loop over all vertices of the entity
      for (size_t local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
        // find the global DoF id of this vertex, therefore
        size_t success = 0;
        size_t failure = 0;
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          assert(ii < std::numeric_limits< int >::max());
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(int(ii));
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            ++success;
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(int(ii));
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            const size_t local_DoF_id = entity_dg_local_key.index();
            const size_t global_DoF_id = range_.space().mapper().mapToGlobal(entity, local_DoF_id);
//            num_local_DoFs_per_global_DoF_[global_DoF_id] += 1;
            // add this global DoF to this vertex
            vertex_to_dof_id_map_[global_vertex_ids_of_current_entity_[local_vertex_id]].insert(global_DoF_id);
//            // evaluate the source and add it to the range
//            const FieldType source_value = local_source_DoF_vector.get(local_DoF_id);
//            range_.vector().add_to_entry(global_DoF_id, source_value);
          } else
            ++failure;
        } // loop over all local DoFs
        assert(success == 1 && failure == (num_vertices - 1)
               && "This must not happen for a polOrder 1 Lagrange space!");
      } // loop over all vertices of the entity
    } // if (!finalized_)
  } // ... apply_local(... EntityType ...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& /*outside_entity*/) /*DS_OVERRIDE DS_FINAL*/
  {
    assert(prepared_);
    if (!finalized_) {
      // some preparations (the griwalker usually calls apply_local(entity) first, so current_entity_ is already
      //                    inside_entity and we do not do any additional work)
      prepare_for_current_(inside_entity);
      const size_t num_vertices = inside_entity.template count< dimDomain >();
      // loop over all vertices of the intersection
      const auto& intersection_geometry = intersection.geometry();
      for (int local_intersection_corner_id = 0;
           local_intersection_corner_id < intersection_geometry.corners();
           ++local_intersection_corner_id) {
        const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
        // now, we need to find the entity's vertex this intersection's corner point equals to, so we
        // loop over all vertices of the entity
        for (size_t local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
          if (Stuff::Common::FloatCmp::eq(global_intersection_corner,
                                          global_vertex_coordinates_of_current_entity_[local_vertex_id]))
            boundary_vertices_.insert(global_vertex_ids_of_current_entity_[local_vertex_id]);
        } // loop over all vertices of the entity
      } // loop over all vertices of the intersection
    } // if (!finalized_)
  } // ... apply_local(... IntersectionType ...)

  virtual void finalize() /*DS_OVERRIDE DS_FINAL*/
  {
    if (!finalized_) {
//      // loop over all DoFs
//      assert(num_local_DoFs_per_global_DoF_.size() == range_.space().mapper().size());
//      for (size_t global_DoF = 0; global_DoF < num_local_DoFs_per_global_DoF_.size(); ++global_DoF) {
//        const FieldType value = range_.vector().get_entry(global_DoF);
//        range_.vector().set_entry(global_DoF, value / num_local_DoFs_per_global_DoF_[global_DoF]);
//      }

      // walk the grid for the second time
      for (auto entity_it = grid_view_.template begin< 0 >(); entity_it != grid_view_.template end< 0 >(); ++entity_it) {
        const auto& entity = *entity_it;
        // get the local functions
        const auto local_source = source_.local_discrete_function(entity);
        const auto& local_source_DoF_vector = local_source.vector();
        // get the local finite elements
        // * for the oswald projection
        typedef typename SpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
        const auto dg_finite_element = range_.space().backend().finiteElement(entity);
        const FiniteElementType cg_finite_element(entity.geometry().type(), 1);
        const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
        const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
        assert(dg_local_coefficients.size() == cg_local_coefficients.size()
               && "Wrong finite element given!");
        // to compute the oswald projection
        // * loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          assert(ii < std::numeric_limits< int >::max());
          const auto& entity_dg_local_key = dg_local_coefficients.localKey(int(ii));
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(int(ii));
          assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
          const int local_vertex_id = entity_cg_local_key.subEntity();
          const size_t local_DoF_id = entity_dg_local_key.index();
          const auto vertexPtr = entity.template subEntity< dimDomain >(local_vertex_id);
          const auto& vertex = *vertexPtr;
          const size_t global_vertex_id = grid_view_.indexSet().index(vertex);
          // if we are on the domain boundary
          if (boundary_vertices_.count(global_vertex_id)) {
            // get global DoF id
            const size_t global_DoF_id = range_.space().mapper().mapToGlobal(entity, local_DoF_id);
            // set the dof to zero (we have dirichlet zero)
            range_.vector().set_entry(global_DoF_id, FieldType(0));
          } else {
            // do the oswald projection
            const size_t num_DoFS_per_vertex = vertex_to_dof_id_map_[global_vertex_id].size();
            // * get the source Dof
            const FieldType source_Dof_value = local_source_DoF_vector.get(local_DoF_id);
            // * and add it to all target DoFs
            for (size_t target_global_DoF_id : vertex_to_dof_id_map_[global_vertex_id])
              range_.vector().add_to_entry(target_global_DoF_id, source_Dof_value / num_DoFS_per_vertex);
          } // if (boundary_vertices.find(global_vertex_id))
        } // loop over all local DoFs
      } // walk the grid for the second time

      // clean up
      cg_fem_ = nullptr;
      dg_fem_ = nullptr;
      global_vertex_ids_of_current_entity_ = std::vector< size_t >();
      global_vertex_coordinates_of_current_entity_ = std::vector< DomainType >();
      vertex_to_dof_id_map_ = std::vector< std::set< size_t > >();
      boundary_vertices_ = std::set< size_t >();
      finalized_ = true;
    } // if (!finalized_)
  } // ... finalize(...)

private:
  void prepare_for_current_(const EntityType& entity)
  {
    const size_t entity_index = grid_view_.indexSet().index(entity);
    if (entity_index != current_entity_) {
      // prepare storage
      const size_t num_vertices = entity.template count< dimDomain >();
      if (global_vertex_ids_of_current_entity_.size() < num_vertices)
        global_vertex_ids_of_current_entity_ = std::vector< size_t >(num_vertices, 0);
      if (global_vertex_coordinates_of_current_entity_.size() < num_vertices)
        global_vertex_coordinates_of_current_entity_ = std::vector< DomainType >(num_vertices, DomainType(0));
      // get the local finite elements
      cg_fem_ = std::unique_ptr< CgFemType >(new CgFemType(entity.geometry().type(), 1));
      dg_fem_ = std::unique_ptr< DgFemType >(new DgFemType(range_.space().backend().finiteElement(entity)));
      const auto& dg_local_coefficients = dg_fem_->localCoefficients();
      const auto& cg_local_coefficients = cg_fem_->localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // loop over all vertices of the entity
      for (size_t local_vertex_id = 0; local_vertex_id < num_vertices; ++local_vertex_id) {
        // get global vertex id and coordinate
        assert(local_vertex_id < std::numeric_limits< int >::max());
        const auto vertexPtr = entity.template subEntity< dimDomain >(int(local_vertex_id));
        const auto& vertex = *vertexPtr;
        global_vertex_ids_of_current_entity_[local_vertex_id] = grid_view_.indexSet().index(vertex);
        global_vertex_coordinates_of_current_entity_[local_vertex_id] = vertex.geometry().center();
      } // loop over all vertices of the entity
    } // if (entity_index != current_entity_)
  } // ... prepare_for_current_(...)

  const GridViewType& grid_view_;
  const SourceType& source_;
  RangeType& range_;
  bool prepared_;
  bool applied_;
  bool finalized_;
  size_t current_entity_;
  std::unique_ptr< const CgFemType > cg_fem_;
  std::unique_ptr< const DgFemType > dg_fem_;
  std::vector< size_t > global_vertex_ids_of_current_entity_;
  std::vector< DomainType > global_vertex_coordinates_of_current_entity_;
  // * a map from a global vertex id to global DoF ids
  //   given a vertex id, one obtains a set of all global DoF ids, which are associated with this vertex
  std::vector< std::set< size_t > > vertex_to_dof_id_map_;
  // * a set to hold the global id of all boundary vertices
  std::vector< size_t > num_local_DoFs_per_global_DoF_;
  std::set< size_t > boundary_vertices_;
}; // class OswaldInterpolationLocalizable
#endif

// forward, to be used in the traits
template <class GridViewImp, class FieldImp = double>
class OswaldInterpolation;


template <class GridViewImp, class FieldImp>
class OswaldInterpolationTraits
{
public:
  typedef OswaldInterpolation<GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class OswaldInterpolationTraits


template <class GridViewImp, class FieldImp>
class OswaldInterpolation : public OperatorInterface<OswaldInterpolationTraits<GridViewImp, FieldImp>>
{
public:
  typedef OswaldInterpolationTraits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  OswaldInterpolation(const GridViewType& grid_view, const bool zero_boundary = true)
    : grid_view_(grid_view)
    , zero_boundary_(zero_boundary)
  {
  }

  template <class SGP, class SV, class RGP, class RV>
  void apply(
      const ConstDiscreteFunction<DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<SGP, 1, FieldType, 1, 1>, SV>&
          source,
      DiscreteFunction<DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<RGP, 1, FieldType, 1, 1>, RV>& range) const
  {
    typedef DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<RGP, 1, FieldType, 1, 1> SpaceType;
    static const unsigned int polOrder = 1;

    // data structures we need
    // * a map from a global vertex id to global DoF ids
    //   given a vertex, one obtains a set of all global DoF ids, which are associated with this vertex
    typedef std::vector<std::set<size_t>> VertexToEntitiesMapType;
    VertexToEntitiesMapType vertex_to_dof_id_map(grid_view_.indexSet().size(dimDomain));
    // * a set to hold the global id off all boundary vertices
    std::set<size_t> boundary_vertices;

    // walk the grid for the first time
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local finite elements
      typedef typename SpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = range.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // loop over all vertices
      std::vector<DomainType> global_vertices(entity.template count<dimDomain>(), DomainType(0));
      std::vector<size_t> global_vertex_ids(global_vertices.size(), 0);
      for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
        // get global vertex id
        const auto vertexPtr               = entity.template subEntity<dimDomain>(int(local_vertex_id));
        const auto& vertex                 = *vertexPtr;
        global_vertex_ids[local_vertex_id] = grid_view_.indexSet().index(vertex);
        global_vertices[local_vertex_id]   = vertex.geometry().center();
        // find the global DoF id to this vertex, therefore
        // loop over all local DoFs
        for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
          assert(ii < std::numeric_limits<int>::max());
          const auto& entity_cg_local_key = cg_local_coefficients.localKey(int(ii));
          if (entity_cg_local_key.subEntity() == local_vertex_id) {
            const auto& entity_dg_local_key = dg_local_coefficients.localKey(int(ii));
            assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
            const size_t local_DOF_id  = entity_dg_local_key.index();
            const size_t global_DOF_id = range.space().mapper().mapToGlobal(entity, local_DOF_id);
            // add this global DoF to this vertex
            vertex_to_dof_id_map[global_vertex_ids[local_vertex_id]].insert(global_DOF_id);
            // there must be one and only one for a polorder 1 lagrange basis
            break;
          }
        } // loop over all local DoFs
      } // loop over all vertices
      if (zero_boundary_) {
        // in order to determine the boundary vertices, we need to
        // loop over all intersections
        const auto intersectionEndIt = grid_view_.iend(entity);
        for (auto intersectionIt = grid_view_.ibegin(entity); intersectionIt != intersectionEndIt; ++intersectionIt) {
          const auto& intersection = *intersectionIt;
          if (intersection.boundary() && !intersection.neighbor()) {
            const auto& intersection_geometry = intersection.geometry();
            for (int local_intersection_corner_id = 0; local_intersection_corner_id < intersection_geometry.corners();
                 ++local_intersection_corner_id) {
              const auto global_intersection_corner = intersection_geometry.corner(local_intersection_corner_id);
              // now, we need to find the entity's vertex this intersection's corner point equals to, so we
              // loop over all vertices of the entity
              for (size_t local_vertex_id = 0; local_vertex_id < global_vertices.size(); ++local_vertex_id) {
                if (Stuff::Common::FloatCmp::eq(global_intersection_corner, global_vertices[local_vertex_id]))
                  boundary_vertices.insert(global_vertex_ids[local_vertex_id]);
              } // loop over all vertices of the entity
            } // if (intersection.boundary() && !intersection.neighbor())
          } // loop over all intersections
        } // loop over all intersections
      }
    } // walk the grid for the first time

    // walk the grid for the second time
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get the local functions
      const auto local_source             = source.local_discrete_function(entity);
      const auto& local_source_DoF_vector = local_source.vector();
      // get the local finite elements
      // * for the oswald projection
      typedef typename SpaceType::Traits::ContinuousFiniteElementType FiniteElementType;
      const auto dg_finite_element = range.space().backend().finiteElement(entity);
      const FiniteElementType cg_finite_element(entity.geometry().type(), polOrder);
      const auto& dg_local_coefficients = dg_finite_element.localCoefficients();
      const auto& cg_local_coefficients = cg_finite_element.localCoefficients();
      assert(dg_local_coefficients.size() == cg_local_coefficients.size() && "Wrong finite element given!");
      // to compute the oswald projection
      // * loop over all local DoFs
      for (size_t ii = 0; ii < dg_local_coefficients.size(); ++ii) {
        assert(ii < std::numeric_limits<int>::max());
        const auto& entity_dg_local_key = dg_local_coefficients.localKey(int(ii));
        const auto& entity_cg_local_key = cg_local_coefficients.localKey(int(ii));
        assert(entity_cg_local_key.codim() == dimDomain && "Wrong finite element given!");
        const int local_vertex_id     = entity_cg_local_key.subEntity();
        const size_t local_DoF_id     = entity_dg_local_key.index();
        const auto vertexPtr          = entity.template subEntity<dimDomain>(local_vertex_id);
        const auto& vertex            = *vertexPtr;
        const size_t global_vertex_id = grid_view_.indexSet().index(vertex);
        // if we are on the domain boundary
        if (zero_boundary_ && boundary_vertices.count(global_vertex_id)) {
          // get global DoF id
          const size_t global_DoF_id = range.space().mapper().mapToGlobal(entity, local_DoF_id);
          // set the dof to zero (we have dirichlet zero)
          range.vector().set_entry(global_DoF_id, FieldType(0));
        } else {
          // do the oswald projection
          const size_t num_DoFS_per_vertex = vertex_to_dof_id_map[global_vertex_id].size();
          // * get the source Dof
          const FieldType source_Dof_value = local_source_DoF_vector.get(local_DoF_id);
          // * and add it to all target DoFs
          for (size_t target_global_DoF_id : vertex_to_dof_id_map[global_vertex_id])
            range.vector().add_to_entry(target_global_DoF_id, source_Dof_value / num_DoFS_per_vertex);
        } // if (boundary_vertices.find(global_vertex_id))
      } // loop over all local DoFs
    } // walk the grid for the second time
  } // ... apply(...)

private:
  const GridViewType& grid_view_;
  const bool zero_boundary_;
}; // class OswaldInterpolation


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_OSWALD_HH
