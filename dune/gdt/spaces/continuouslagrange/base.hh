// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH

#include <dune/common/dynvector.hh>
#include <dune/common/version.hh>
#include <dune/common/deprecated.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,3,9) //EXADUNE
# include <dune/geometry/referenceelements.hh>
#else
# include <dune/geometry/genericreferenceelements.hh>
#endif

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {


// forward, to allow for specialization
template< class ImpTraits, size_t domainDim, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols = 1 >
class
    DUNE_DEPRECATED_MSG("Include <dune/gdt/spaces/cg.hh> and use CGInterface instead (21.11.2014)!")
      ContinuousLagrangeBase
{
  static_assert(AlwaysFalse< ImpTraits >::value, "Untested for these dimensions!");
};


template< class ImpTraits, size_t domainDim, class RangeFieldImp, size_t rangeDim >
class
    DUNE_DEPRECATED_MSG("Include <dune/gdt/spaces/cg.hh> and use CGInterface instead (21.11.2014)!")
      ContinuousLagrangeBase< ImpTraits, domainDim, RangeFieldImp, rangeDim, 1 >
  : public SpaceInterface< ImpTraits, domainDim, rangeDim, 1 >
{
  typedef SpaceInterface< ImpTraits, domainDim, rangeDim, 1 > BaseType;
  typedef ContinuousLagrangeBase< ImpTraits, domainDim, RangeFieldImp, rangeDim, 1 > ThisType;

  static constexpr RangeFieldImp compare_tolerance_ = 1e-13;
public:
  typedef ImpTraits Traits;

  using BaseType::polOrder;

  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::DomainType;

  typedef typename Traits::RangeFieldType RangeFieldType;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;

  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  using typename BaseType::BoundaryInfoType;
  using typename BaseType::PatternType;

  virtual ~ContinuousLagrangeBase() {}

  using BaseType::compute_pattern;

  template< class G, class S, size_t d, size_t r, size_t rC >
  PatternType compute_pattern(const GridView< G >& local_grid_view, const SpaceInterface< S, d, r, rC >& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(local_grid_view, ansatz_space);
  }

  virtual std::vector< DomainType > lagrange_points(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1) DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view().indexSet().contains(entity));
    // get the basis and reference element
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector< RangeType > tmp_basis_values(basis.size(), RangeType(0));
    const auto& reference_element = ReferenceElements< DomainFieldType, dimDomain >::general(entity.type());
    const int num_vertices = reference_element.size(dimDomain);
    assert(num_vertices >= 0);
    assert(size_t(num_vertices) == basis.size() && "This should not happen with polOrder 1!");
    // prepare return vector
    std::vector< DomainType > local_vertices(num_vertices, DomainType(0));
    // loop over all vertices
    for (int ii = 0; ii < num_vertices; ++ii) {
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
  } // ... lagrange_points(...)

  virtual std::set< size_t > local_dirichlet_DoFs(const EntityType& entity,
                                                  const BoundaryInfoType& boundaryInfo) const
  {
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1) DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    // check
    assert(this->grid_view().indexSet().contains(entity));
    // prepare
    std::set< size_t > localDirichletDofs;
    std::vector< DomainType > dirichlet_vertices;
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = this->grid_view().iend(entity);
    for (auto intersection_it = this->grid_view().ibegin(entity);
         intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundaryInfo.dirichlet(intersection) || (!intersection.neighbor() && !intersection.boundary())) {
        // and get the vertices of the intersection
        const auto geometry = intersection.geometry();
        for (int cc = 0; cc < geometry.corners(); ++cc)
          dirichlet_vertices.emplace_back(entity.geometry().local(geometry.corner(cc)));
      } // only work on dirichlet ones
    } // loop over all intersections
    // find the corresponding basis functions
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector< RangeType > tmp_basis_values(basis.size(), RangeType(0));
    for (size_t cc = 0; cc < dirichlet_vertices.size(); ++cc) {
      // find the basis function that evaluates to one here (has to be only one!)
      basis.evaluate(dirichlet_vertices[cc], tmp_basis_values);
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
  } // ... local_dirichlet_DoFs(...)

  using BaseType::local_constraints;

  template< class S, size_t d, size_t r, size_t rC, class ConstraintsType >
  void local_constraints(const SpaceInterface< S, d, r, rC >& /*other*/,
                         const EntityType& /*entity*/,
                         ConstraintsType& /*ret*/) const
  {
    static_assert(AlwaysFalse< S >::value, "Not implemented for these constraints!");
  }

  template< class S, size_t d, size_t r, size_t rC >
  void local_constraints(const SpaceInterface< S, d, r, rC >& other,
                         const EntityType& entity,
                         DirichletConstraints< IntersectionType >& ret) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1) DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view().indexSet().contains(entity));
    const std::set< size_t > localDirichletDofs = this->local_dirichlet_DoFs(entity, ret.boundary_info());
    const size_t numRows = localDirichletDofs.size();
    Dune::DynamicVector< size_t > tmpMappedRows;
    Dune::DynamicVector< size_t > tmpMappedCols;
    if (numRows > 0) {
      const size_t numCols = this->mapper().numDofs(entity);
      ret.set_size(numRows, numCols);
      this->mapper().globalIndices(entity, tmpMappedRows);
      other.mapper().globalIndices(entity, tmpMappedCols);
      size_t localRow = 0;
      for (const size_t& localDirichletDofIndex : localDirichletDofs) {
        ret.global_row(localRow) = tmpMappedRows[localDirichletDofIndex];
        for (size_t jj = 0; jj < ret.cols(); ++jj) {
          ret.global_col(jj) = tmpMappedCols[jj];
          if (tmpMappedCols[jj] == tmpMappedRows[localDirichletDofIndex])
            ret.value(localRow, jj) = ret.set_row() ? 1 : 0;
          else
            ret.value(localRow, jj) = 0;
        }
        ++localRow;
      }
    } else {
      ret.set_size(0, 0);
    }
  } // ... local_constraints(..., Constraints::Dirichlet< ... > ...)
}; // class ContinuousLagrangeBase


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH
