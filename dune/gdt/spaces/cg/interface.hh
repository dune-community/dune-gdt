// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CG_INTERFACE_HH
#define DUNE_GDT_SPACES_CG_INTERFACE_HH

#include <algorithm>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/typetraits.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) //EXADUNE
# include <dune/geometry/referenceelements.hh>
#else
# include <dune/geometry/genericreferenceelements.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/float_cmp.hh>
#include <dune/stuff/common/type_utils.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {


static constexpr ChooseSpaceBackend default_cg_backend = default_space_backend;


template< class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1 >
class CGInterface
  : public SpaceInterface< ImpTraits, domainDim, rangeDim, rangeDimCols >
{
  typedef SpaceInterface< ImpTraits, domainDim, rangeDim, rangeDimCols > BaseType;
  typedef CGInterface< ImpTraits, domainDim, rangeDim, rangeDimCols > ThisType;
public:
  typedef ImpTraits Traits;

  using BaseType::polOrder;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;

  using typename BaseType::RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;

  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::PatternType;
private:
  static const constexpr RangeFieldType compare_tolerance_ = 1e-13;
public:

  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/
  std::vector< DomainType > lagrange_points(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().lagrange_points(entity));
    return this->as_imp().lagrange_points(entity);
  } // ... lagrange_points(...)

  std::set< size_t > local_dirichlet_DoFs(const EntityType& entity,
                                          const BoundaryInfoType& boundaryInfo) const
  {
    CHECK_CRTP(this->as_imp().local_dirichlet_DoFs(entity, boundaryInfo));
    return this->as_imp().local_dirichlet_DoFs(entity, boundaryInfo);
  } // ... local_dirichlet_DoFs(...)
  /** @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience.''
   * @{
   **/

  /**
   * \todo Use our FloatCmp with the correct comparison method gainst 0 and 1!
   */
  std::vector< DomainType > lagrange_points_order_1(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Does not work for higher polynomial orders!");
    if (dimRange != 1) DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view().indexSet().contains(entity));
    // get the basis and reference element
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector< RangeType > tmp_basis_values(basis.size(), RangeType(0));
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity.type());
    const auto num_vertices = reference_element.size(dimDomain);
    assert(num_vertices >= 0);
    assert(boost::numeric_cast< size_t >(num_vertices) == basis.size() && "This should not happen with polOrder 1!");
    // prepare return vector
    std::vector< DomainType > local_vertices(num_vertices, DomainType(0));
    // loop over all vertices
    for (auto ii : DSC::valueRange(num_vertices)) {
      // get the local coordinate of the iith vertex
      const auto local_vertex = reference_element.position(ii, dimDomain);
      // evaluate the basefunctionset
      basis.evaluate(local_vertex, tmp_basis_values);
      // find the basis function that evaluates to one here (has to be only one!)
      size_t ones = 0;
      size_t zeros = 0;
      size_t failures = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj) {
        if (std::abs((tmp_basis_values)[jj][0] - RangeFieldType(1)) < compare_tolerance_) {
          local_vertices[jj] = local_vertex;
          ++ones;
        } else if (std::abs((tmp_basis_values)[jj][0]) < compare_tolerance_)
          ++zeros;
        else
          ++failures;
      }
      assert(ones == 1 && zeros == (basis.size() - 1) && failures == 0 && "This must not happen for polOrder 1!");
    }
    return local_vertices;
  } // ... lagrange_points_order_1(...)

  /**
   * \attention The current search strategy is not rubust with respect to the grid (see note)!
   * \note      We only look in the intersections of entity and its direct neighbors for dirichlet vertices and might
   *            miss some (see also https://github.com/pymor/dune-gdt/issues/51)!
   */
  std::set< size_t > local_dirichlet_DoFs_order_1(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1) DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    // check
    assert(this->grid_view().indexSet().contains(entity));
    // prepare
    std::set< size_t > localDirichletDofs;
    std::vector< DomainType > global_dirichlet_vertices;
    // find all dirichlet vertices of this entity and its neighbors (in global coordinates)
    add_dirichlet_vertices(boundaryInfo, entity, /*recursion_level = */ 1, global_dirichlet_vertices);
    // keep those local to this entity
    std::vector< DomainType > dirichlet_vertices;
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity.type());
    for (const auto& global_vertex : global_dirichlet_vertices) {
      auto local_vertex = entity.geometry().local(global_vertex);
      if (reference_element.checkInside(local_vertex) && std::find(dirichlet_vertices.begin(),
                                                                   dirichlet_vertices.end(),
                                                                   local_vertex) == dirichlet_vertices.end())
        dirichlet_vertices.push_back(local_vertex);
    }
    // find the corresponding basis functions
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector< RangeType > tmp_basis_values(basis.size(), RangeType(0));
    for (const auto& dirichlet_vertex : dirichlet_vertices) {
      // find the basis function that evaluates to one here (has to be only one!)
      basis.evaluate(dirichlet_vertex, tmp_basis_values);
      size_t ones = 0;
      size_t zeros = 0;
      size_t failures = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj) {
        if (std::abs(tmp_basis_values[jj][0] - RangeFieldType(1)) < compare_tolerance_) {
          localDirichletDofs.insert(jj);
          ++ones;
        } else if (std::abs(tmp_basis_values[jj][0]) < compare_tolerance_)
          ++zeros;
        else
          ++failures;
      }
      assert(ones == 1 && zeros == (basis.size() - 1) && failures == 0 && "This must not happen for polOrder 1!");
    }
    return localDirichletDofs;
  } // ... local_dirichlet_DoFs_order_1(...)

  using BaseType::compute_pattern;

  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(local_grid_view, ansatz_space);
  }

  using BaseType::local_constraints;

  template< class S, size_t d, size_t r, size_t rC, class ConstraintsType >
  void local_constraints(const SpaceInterface< S, d, r, rC >& /*other*/,
                         const EntityType& /*entity*/,
                         ConstraintsType& /*ret*/) const
  {
    static_assert(AlwaysFalse< S >::value, "Not implemented for these constraints!");
  }

  template< class S, size_t d, size_t r, size_t rC >
  void local_constraints(const SpaceInterface< S, d, r, rC >& /*other*/,
                         const EntityType& entity,
                         DirichletConstraints< IntersectionType >& ret) const
  {
    const auto local_DoFs = this->local_dirichlet_DoFs(entity, ret.boundary_info());
    if (local_DoFs.size() > 0) {
      const auto global_indices = this->mapper().globalIndices(entity);
      for (const auto& local_DoF : local_DoFs) {
        ret.insert(global_indices[local_DoF]);
      }
    }
  } // ... local_constraints(..., Constraints::Dirichlet< ... > ...)
  /** @} */

private:
  void add_dirichlet_vertices(const BoundaryInfoType& boundaryInfo,
                              const EntityType& entity,
                              const ssize_t recursion_level,
                              std::vector< DomainType >& dirichlet_vertices) const
  {
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = this->grid_view().iend(entity);
    for (auto intersection_it = this->grid_view().ibegin(entity);
         intersection_it != intersection_it_end;
         ++intersection_it) {
      const auto& intersection = *intersection_it;
      // only work on dirichlet intersections + process boundaries for parallel runs
      if (boundaryInfo.dirichlet(intersection) || (!intersection.neighbor() && !intersection.boundary())) {
        // and get the vertices of the intersection
        const auto geometry = intersection.geometry();
        for (auto cc : DSC::valueRange(geometry.corners()))
          dirichlet_vertices.emplace_back(geometry.corner(cc));
      } // only work on dirichlet intersections + process boundaries for parallel runs
      if (recursion_level > 0) {
        // also call myself on all neighbors
        if (intersection.neighbor()) {
          const auto neighbor_ptr = intersection.outside();
          add_dirichlet_vertices(boundaryInfo, *neighbor_ptr, recursion_level - 1, dirichlet_vertices);
        }
      } // if (level > 0)
    } // loop over all intersections
  } // ... add_dirichlet_vertices(...)
}; // class CGInterface


} // namespace Spaces
namespace internal {


template< class S >
struct is_cg_space_helper
{
  DSC_has_typedef_initialize_once(Traits)
  DSC_has_static_member_initialize_once(dimDomain)
  DSC_has_static_member_initialize_once(dimRange)
  DSC_has_static_member_initialize_once(dimRangeCols)

  static const bool is_candidate = DSC_has_typedef(Traits)< S >::value
                                   && DSC_has_static_member(dimDomain)< S >::value
                                   && DSC_has_static_member(dimRange)< S >::value
                                   && DSC_has_static_member(dimRangeCols)< S >::value;
}; // class is_cg_space_helper


} // namespace internal


template< class S, bool candidate = internal::is_cg_space_helper< S >::is_candidate >
struct is_cg_space
  : public std::is_base_of< Spaces::CGInterface< typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols >
                          , S >
{};


template< class S >
struct is_cg_space< S, false >
  : public std::false_type
{};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_INTERFACE_HH
