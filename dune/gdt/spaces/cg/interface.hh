// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_CG_INTERFACE_HH
#define DUNE_GDT_SPACES_CG_INTERFACE_HH

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/common/ranges.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {


static constexpr ChooseSpaceBackend default_cg_backend = default_space_backend;


template <class ImpTraits, size_t domainDim, size_t rangeDim, size_t rangeDimCols = 1>
class CgSpaceInterface : public SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols>
{
  typedef SpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols> BaseType;
  typedef CgSpaceInterface<ImpTraits, domainDim, rangeDim, rangeDimCols> ThisType;

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
  typedef Dune::XT::Common::FieldVector<DomainFieldType, dimDomain> StuffDomainType;
  static const constexpr RangeFieldType compare_tolerance_ = 1e-13;

public:
  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/
  std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().lagrange_points(entity));
    return this->as_imp().lagrange_points(entity);
  } // ... lagrange_points(...)

  std::set<size_t> local_dirichlet_DoFs(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    CHECK_CRTP(this->as_imp().local_dirichlet_DoFs(entity, boundaryInfo));
    return this->as_imp().local_dirichlet_DoFs(entity, boundaryInfo);
  } // ... local_dirichlet_DoFs(...)
  /** @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience.''
   * @{
   **/
  std::vector<DomainType> lagrange_points_order_1(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1)
      DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view().indexSet().contains(entity));
    // get the basis and reference element
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector<RangeType> tmp_basis_values(basis.size(), RangeType(0));
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity.type());
    const auto num_vertices       = reference_element.size(dimDomain);
    assert(num_vertices >= 0);
    assert(boost::numeric_cast<size_t>(num_vertices) == basis.size() && "This should not happen with polOrder 1!");
    // prepare return vector
    std::vector<DomainType> local_vertices(num_vertices, DomainType(0));
    // loop over all vertices
    for (auto ii : Dune::XT::Common::value_range(num_vertices)) {
      // get the local coordinate of the iith vertex
      const auto local_vertex = reference_element.position(ii, dimDomain);
      // evaluate the basefunctionset
      basis.evaluate(local_vertex, tmp_basis_values);
      // find the basis function that evaluates to one here (has to be only one!)
      size_t ones     = 0;
      size_t zeros    = 0;
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

  std::set<size_t> local_dirichlet_DoFs_order_1(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    static const XT::Grid::DirichletBoundary dirichlet{};
    if (dimRange != 1)
      DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    // check
    assert(this->grid_view().indexSet().contains(entity));
    if (!entity.hasBoundaryIntersections())
      return std::set<size_t>();
    // prepare
    std::set<size_t> localDirichletDofs;
    std::vector<DomainType> dirichlet_vertices;
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = this->grid_view().iend(entity);
    for (auto intersection_it = this->grid_view().ibegin(entity); intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundaryInfo.type(intersection) == dirichlet || (!intersection.neighbor() && !intersection.boundary())) {
        // and get the vertices of the intersection
        const auto geometry = intersection.geometry();
        for (auto cc : Dune::XT::Common::value_range(geometry.corners()))
          dirichlet_vertices.emplace_back(entity.geometry().local(geometry.corner(cc)));
      } // only work on dirichlet ones
    } // loop over all intersections
    // find the corresponding basis functions
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector<RangeType> tmp_basis_values(basis.size(), RangeType(0));
    for (size_t cc = 0; cc < dirichlet_vertices.size(); ++cc) {
      // find the basis function that evaluates to one here (has to be only one!)
      basis.evaluate(dirichlet_vertices[cc], tmp_basis_values);
      size_t ones     = 0;
      size_t zeros    = 0;
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
      // asserts valid for polorder 1 only
      assert(ones == 1);
      assert(zeros == (basis.size() - 1));
      assert(failures == 0);
    }
    return localDirichletDofs;
  } // ... local_dirichlet_DoFs_order_1(...)


  std::set<size_t> local_dirichlet_DoFs_simplicial_lagrange_elements(const EntityType& entity,
                                                                     const BoundaryInfoType& boundaryInfo) const
  {
    if (!entity.type().isSimplex())
      DUNE_THROW(NotImplemented, "Only implemented for simplex elements!");
    // check
    assert(this->grid_view().indexSet().contains(entity));
    // prepare
    std::set<size_t> localDirichletDofs;
    std::vector<DomainType> dirichlet_vertices;
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = this->grid_view().iend(entity);
    for (auto intersection_it = this->grid_view().ibegin(entity); intersection_it != intersection_it_end;
         ++intersection_it) {
      const auto& intersection = *intersection_it;
      std::vector<StuffDomainType> dirichlet_vertices_intersection;
      // only work on dirichlet ones, actual dirichlet intersections + process boundaries for parallel runs
      if (boundaryInfo.dirichlet(intersection) || (!intersection.neighbor() && !intersection.boundary())) {
        // and get the vertices of the intersection
        const auto geometry = intersection.geometry();
        for (auto cc : Dune::XT::Common::value_range(geometry.corners())) {
          dirichlet_vertices_intersection.emplace_back(entity.geometry().local(geometry.corner(cc)));
          dirichlet_vertices.emplace_back(entity.geometry().local(geometry.corner(cc)));
        }
      } // only work on dirichlet ones
      // for higher polynomial orders, add points associated to DoFs within the intersection
      // calculate by using lagrange grid {x = \sum_{j=0}^d \lambda_j a_j | \sum_j lambda_j = 1}, where a_j are the
      // vertices of the entity and \lambda_j \in {\frac{m}{polOrder} | m = 0, ... , polOrder}
      std::vector<double> possible_coefficients(polOrder < 1 ? 0 : polOrder - 1);
      for (int m                 = 0; m < polOrder; ++m)
        possible_coefficients[m] = m / polOrder;
      std::set<std::vector<double>> possible_coefficient_vectors;
      possible_convex_combination_coefficients(
          possible_coefficient_vectors, possible_coefficients, dirichlet_vertices_intersection.size());
      for (const auto& coefficient_vector : possible_coefficient_vectors)
        dirichlet_vertices.emplace_back(std::inner_product(dirichlet_vertices_intersection.begin(),
                                                           dirichlet_vertices_intersection.end(),
                                                           coefficient_vector.begin(),
                                                           StuffDomainType(0)));
    } // loop over all intersections

    // find the corresponding basis functions
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector<RangeType> tmp_basis_values(basis.size(), RangeType(0));
    for (size_t cc = 0; cc < dirichlet_vertices.size(); ++cc) {
      // find the basis function that evaluates to one here (has to be only one per range dimension!)
      basis.evaluate(dirichlet_vertices[cc], tmp_basis_values);
      size_t ones     = 0;
      size_t zeros    = 0;
      size_t failures = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj) {
        for (size_t rr = 0; rr < dimRange; ++rr) {
          if (std::abs(tmp_basis_values[jj][rr] - RangeFieldType(1)) < compare_tolerance_) {
            localDirichletDofs.insert(jj);
            ++ones;
          } else if (std::abs(tmp_basis_values[jj][rr]) < compare_tolerance_)
            ++zeros;
          else
            ++failures;
        }
      }
      assert(ones == dimRange && zeros == ((basis.size() - 1) * dimRange) && failures == 0 && "This must not happen!");
    }
    return localDirichletDofs;
  } // ... local_dirichlet_DoFs_simplicial_lagrange_elements(...)

  using BaseType::compute_pattern;

  template <class G, class S, size_t d, size_t r, size_t rC>
  PatternType compute_pattern(const GridView<G>& local_grid_view, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(local_grid_view, ansatz_space);
  }

  using BaseType::local_constraints;

  template <class S, size_t d, size_t r, size_t rC, class ConstraintsType>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*other*/, const EntityType& /*entity*/,
                         ConstraintsType& /*ret*/) const
  {
    static_assert(AlwaysFalse<S>::value, "Not implemented for these constraints!");
  }

  template <class S, size_t d, size_t r, size_t rC>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*other*/, const EntityType& entity,
                         DirichletConstraints<IntersectionType>& ret) const
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
  void possible_convex_combination_coefficients(std::set<std::vector<double>>& vectors_in,
                                                const std::vector<double>& possible_coefficients,
                                                const size_t final_size) const
  {
    if (vectors_in.empty()) {
      if (final_size != 0)
        for (const auto& coeff : possible_coefficients)
          vectors_in.insert(std::vector<double>(1, coeff));
    } else if (vectors_in.begin()->size() != final_size) {
      std::set<std::vector<double>> vectors_out;
      for (auto& vec : vectors_in) {
        for (const auto& coeff : possible_coefficients) {
          if ((vec.size() != final_size - 1
               && Dune::XT::Common::FloatCmp::le(std::accumulate(vec.begin(), vec.end(), 0.0) + coeff, 1.0))
              || (vec.size() == final_size - 1
                  && Dune::XT::Common::FloatCmp::eq(std::accumulate(vec.begin(), vec.end(), 0.0) + coeff, 1.0))) {
            std::vector<double> vec_copy = vec;
            vec_copy.push_back(coeff);
            vectors_out.insert(vec_copy);
          }
        }
      }
      vectors_in = vectors_out;
      possible_convex_combination_coefficients(vectors_in, possible_coefficients, final_size);
    } // if (...)
  } // ... possible_convex_combination_coefficients(...)
}; // class CgSpaceInterface


namespace internal {


template <class S>
struct is_cg_space_helper
{
  DXTC_has_typedef_initialize_once(Traits) DXTC_has_static_member_initialize_once(dimDomain)
      DXTC_has_static_member_initialize_once(dimRange) DXTC_has_static_member_initialize_once(dimRangeCols)

          static const
      bool is_candidate = DXTC_has_typedef(Traits)<S>::value && DXTC_has_static_member(dimDomain)<S>::value
                          && DXTC_has_static_member(dimRange)<S>::value
                          && DXTC_has_static_member(dimRangeCols)<S>::value;
}; // class is_cg_space_helper


} // namespace internal


template <class S, bool candidate = internal::is_cg_space_helper<S>::is_candidate>
struct is_cg_space
    : public std::is_base_of<CgSpaceInterface<typename S::Traits, S::dimDomain, S::dimRange, S::dimRangeCols>, S>
{
};


template <class S>
struct is_cg_space<S, false> : public std::false_type
{
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_INTERFACE_HH
