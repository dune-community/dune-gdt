// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_RT_INTERFACE_HH
#define DUNE_GDT_SPACES_RT_INTERFACE_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/xt/common/timedlogging.hh>

#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1>
class RtSpaceInterface : public SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols>
{
  typedef SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols> BaseType;
  typedef RtSpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols> ThisType;

public:
  typedef ImpTraits Traits;

  using BaseType::polOrder;
  using BaseType::dimDomain;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::EntityType;
  using typename BaseType::BaseFunctionSetType;
  using typename BaseType::PatternType;

  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/

  /**
   *  \brief  Computes a vector 'indices' of length entity.count< 1 >(), where 'indices[intersection.indexInInside()]'
   *          is the index of the basis function (aka the local DoF index) corresponding to the intersection.
   */
  std::vector<size_t> local_DoF_indices(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().local_DoF_indices(entity));
    return this->as_imp().local_DoF_indices(entity);
  }

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience.''
   * @{
   **/
  std::vector<size_t> local_DoF_indices_2dsimplex_order0(const EntityType& entity) const
  {
    static_assert(dimDomain == 2, "Not implemented!");
    static_assert(polOrder == 0, "Not implemented!");
    // prepare
    const auto num_intersections = entity.subEntities(1);
    std::vector<size_t> local_DoF_index_of_vertex(num_intersections, std::numeric_limits<size_t>::infinity());
    std::vector<size_t> local_DoF_index_of_intersection(num_intersections, std::numeric_limits<size_t>::infinity());
    typedef typename BaseFunctionSetType::DomainType DomainType;
    std::vector<DomainType> vertices(num_intersections, DomainType(0));
    std::vector<bool> lies_on_intersection(num_intersections, false);
    DomainType corner(0);
    typedef typename BaseFunctionSetType::RangeType RangeType;
    const RangeType one(1);
    std::vector<RangeType> basis_values(num_intersections, one);
    const auto basis = this->base_function_set(entity);
    assert(basis.size() == num_intersections);
    const auto geometry = entity.geometry();
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(geometry.type());
    // find the basis function index that corresponds to each vertex of the entity
    // (find the basis function that evaluates to zero at the vertex, and nonzero at the other ones)
    // therefore we walk the vertices
    assert(num_intersections == entity.subEntities(dimDomain));
    for (size_t vv = 0; vv < num_intersections; ++vv) {
      const auto vertex = entity.template subEntity<dimDomain>(boost::numeric_cast<int>(vv));
      // get the vertex coordinates
      vertices[vv] = vertex.geometry().center();
      const auto& vertex_entity = reference_element.position(boost::numeric_cast<int>(vv), dimDomain);
      // evalaute the basis
      basis.evaluate(vertex_entity, basis_values);
      // and find the basis that evaluates zero here
      size_t zeros = 0;
      size_t nonzeros = 0;
      for (size_t ii = 0; ii < num_intersections; ++ii) {
        // we would like to check against 0, but there is a bug in dune-commons FloatCmp
        if (XT::Common::FloatCmp::eq(basis_values[ii] + one, one)) {
          // this is a candidate for the basis function we are looking for
          local_DoF_index_of_vertex[vv] = ii;
          ++zeros;
        } else
          ++nonzeros;
      }
      // make sure there was only one candidate
      if (zeros != 1 || nonzeros != (num_intersections - 1))
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "This must not happen for RTN0 in 2d!\n"
                       << "  zeros    = "
                       << zeros
                       << "\n"
                       << "  nonzeros = "
                       << nonzeros);
    } // walk the vertices
    // so from here on we have the local DoF index that corresponds to each vertex vv in local_DoF_index_of_vertex[vv]
    // now we need to find the intersection opposite to this vertex
    // therefore we walk the intersections
    size_t intersection_counter = 0;
    const auto intersection_it_end = this->grid_layer().iend(entity);
    for (auto intersection_it = this->grid_layer().ibegin(entity); intersection_it != intersection_it_end;
         ++intersection_it) {
      const auto& intersection = *intersection_it;
      const auto intersection_geometry = intersection.geometry();
      const size_t local_intersection_index = intersection.indexInInside();
      // make sure this index has not been already taken by another intersection
      assert(local_DoF_index_of_intersection[local_intersection_index] == std::numeric_limits<size_t>::infinity());
      // walk the corners of the intersection
      for (int cc = 0; cc < intersection_geometry.corners(); ++cc) {
        corner = intersection_geometry.corner(boost::numeric_cast<int>(cc));
        // check which vertices lie on the intersection
        for (size_t vv = 0; vv < num_intersections; ++vv)
          if (XT::Common::FloatCmp::eq(vertices[vv], corner))
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
        DUNE_THROW(XT::Common::Exceptions::internal_error,
                   "This must not happen for RTN0 in 2d!\n"
                       << "  found  = "
                       << found
                       << "\n"
                       << "  missed = "
                       << missed);
      ++intersection_counter;
    } // walk the intersection
    assert(intersection_counter == num_intersections);
    return local_DoF_index_of_intersection;
  } // ... local_DoF_indices_2dsimplex_order0(...)

  using BaseType::compute_pattern;

  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    Dune::XT::Common::TimedLogger().get("gdt.spaces.rt.compute_pattern").warn() << "Returning largest possible pattern!"
                                                                                << std::endl;
    return BaseType::compute_face_and_volume_pattern(grd_layr, ansatz_space);
  }

  using BaseType::local_constraints;

  template <class S, size_t d, size_t r, size_t rC, class C, class R>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*ansatz_space*/,
                         const EntityType& /*entity*/,
                         ConstraintsInterface<C>& /*ret*/) const
  {
    DUNE_THROW(NotImplemented, "RT spaces do not implement constraints!");
  }

  static constexpr bool associates_data_with(int codim)
  {
    return codim < dimDomain;
  }
  /** @} */
}; // class RtSpaceInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_RT_INTERFACE_HH
