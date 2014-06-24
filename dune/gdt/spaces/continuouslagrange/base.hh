// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH
#define DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH

#include <type_traits>

#include <dune/common/typetraits.hh>
#include <dune/common/dynvector.hh>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/geometry/genericreferenceelements.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/exceptions.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {
namespace Spaces {


// forward, to allow for specialization
template <class ImpTraits, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1>
class ContinuousLagrangeBase
{
  static_assert(AlwaysFalse<ImpTraits>::value, "Untested for these dimensions!");
};


template <class ImpTraits, int domainDim, class RangeFieldImp, int rangeDim>
class ContinuousLagrangeBase<ImpTraits, domainDim, RangeFieldImp, rangeDim, 1> : public SpaceInterface<ImpTraits>
{
  typedef SpaceInterface<ImpTraits> BaseType;
  typedef ContinuousLagrangeBase<ImpTraits, domainDim, RangeFieldImp, rangeDim, 1> ThisType;

  static const constexpr RangeFieldImp compare_tolerance_ = 1e-13;

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

  virtual ~ContinuousLagrangeBase()
  {
  }

  using BaseType::compute_pattern;

  template <class G, class S>
  PatternType compute_pattern(const GridView<G>& local_grid_view, const SpaceInterface<S>& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(local_grid_view, ansatz_space);
  }

  virtual std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1)
      DUNE_THROW_COLORFULLY(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view()->indexSet().contains(entity));
    // get the basis and reference element
    const auto basis              = this->base_function_set(entity);
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity.type());
    const int num_vertices = reference_element.size(dimDomain);
    assert(num_vertices >= 0);
    assert(size_t(num_vertices) == basis.size() && "This should not happen with polOrder 1!");
    // prepare return vector
    std::vector<DomainType> local_vertices(num_vertices, DomainType(0));
    if (this->tmp_basis_values_.size() < basis.size())
      this->tmp_basis_values_.resize(basis.size());
    // loop over all vertices
    for (int ii = 0; ii < num_vertices; ++ii) {
      // get the local coordinate of the iith vertex
      const auto local_vertex = reference_element.position(ii, dimDomain);
      // evaluate the basefunctionset
      basis.evaluate(local_vertex, this->tmp_basis_values_);
      // find the basis function that evaluates to one here (has to be only one!)
      size_t ones     = 0;
      size_t zeros    = 0;
      size_t failures = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj) {
        if (std::abs(this->tmp_basis_values_[jj][0] - RangeFieldType(1)) < compare_tolerance_) {
          local_vertices[jj] = local_vertex;
          ++ones;
        } else if (std::abs(this->tmp_basis_values_[jj][0]) < compare_tolerance_)
          ++zeros;
        else
          ++failures;
      }
      assert(ones == 1 && zeros == (basis.size() - 1) && failures == 0 && "This must not happen for polOrder 1!");
    }
    return local_vertices;
  } // ... lagrange_points(...)

  virtual std::set<size_t> local_dirichlet_DoFs(const EntityType& entity, const BoundaryInfoType& boundaryInfo) const
  {
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1)
      DUNE_THROW_COLORFULLY(NotImplemented, "Does not work for higher dimensions");
    // check
    assert(this->grid_view()->indexSet().contains(entity));
    // prepare
    std::set<size_t> localDirichletDofs;
    std::vector<DomainType> dirichlet_vertices;
    // get all dirichlet vertices of this entity, therefore
    // * loop over all intersections
    const auto intersection_it_end = this->grid_view()->iend(entity);
    for (auto intersection_it = this->grid_view()->ibegin(entity); intersection_it != intersection_it_end;
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
    if (this->tmp_basis_values_.size() < basis.size())
      this->tmp_basis_values_.resize(basis.size());
    for (size_t cc = 0; cc < dirichlet_vertices.size(); ++cc) {
      // find the basis function that evaluates to one here (has to be only one!)
      basis.evaluate(dirichlet_vertices[cc], this->tmp_basis_values_);
      size_t ones     = 0;
      size_t zeros    = 0;
      size_t failures = 0;
      for (size_t jj = 0; jj < basis.size(); ++jj) {
        if (std::abs(this->tmp_basis_values_[jj][0] - RangeFieldType(1)) < compare_tolerance_) {
          localDirichletDofs.insert(jj);
          ++ones;
        } else if (std::abs(this->tmp_basis_values_[jj][0]) < compare_tolerance_)
          ++zeros;
        else
          ++failures;
      }
      assert(ones == 1 && zeros == (basis.size() - 1) && failures == 0 && "This must not happen for polOrder 1!");
    }
    return localDirichletDofs;
  } // ... local_dirichlet_DoFs(...)

private:
  template <class C, bool set_row>
  struct DirichletConstraints;
  template <class C>
  struct DirichletConstraints<C, true>
  {
    static RangeFieldType value()
    {
      return RangeFieldType(1);
    }
  };
  template <class C>
  struct DirichletConstraints<C, false>
  {
    static RangeFieldType value()
    {
      return RangeFieldType(0);
    }
  };

  template <class T, bool set_row>
  void compute_local_constraints(const SpaceInterface<T>& other, const EntityType& entity,
                                 Constraints::Dirichlet<IntersectionType, RangeFieldType, set_row>& ret) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1)
      DUNE_THROW_COLORFULLY(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_view()->indexSet().contains(entity));
    typedef DirichletConstraints<Constraints::Dirichlet<IntersectionType, RangeFieldType, set_row>, set_row> SetRow;
    const std::set<size_t> localDirichletDofs = this->local_dirichlet_DoFs(entity, ret.gridBoundary());
    const size_t numRows = localDirichletDofs.size();
    if (numRows > 0) {
      const size_t numCols = this->mapper().numDofs(entity);
      ret.setSize(numRows, numCols);
      this->mapper().globalIndices(entity, tmpMappedRows_);
      other.mapper().globalIndices(entity, tmpMappedCols_);
      size_t localRow = 0;
      const RangeFieldType zero(0);
      for (auto localDirichletDofIt = localDirichletDofs.begin(); localDirichletDofIt != localDirichletDofs.end();
           ++localDirichletDofIt) {
        const size_t& localDirichletDofIndex = *localDirichletDofIt;
        ret.globalRow(localRow) = tmpMappedRows_[localDirichletDofIndex];
        for (size_t jj = 0; jj < ret.cols(); ++jj) {
          ret.globalCol(jj) = tmpMappedCols_[jj];
          if (tmpMappedCols_[jj] == tmpMappedRows_[localDirichletDofIndex])
            ret.value(localRow, jj) = SetRow::value();
          else
            ret.value(localRow, jj) = zero;
        }
        ++localRow;
      }
    } else {
      ret.setSize(0, 0);
    }
  } // ... compute_local_constraints(..., Dirichlet< ..., true >)

public:
  template <bool set>
  void local_constraints(const EntityType& entity,
                         Constraints::Dirichlet<IntersectionType, RangeFieldType, set>& ret) const
  {
    local_constraints(*this, entity, ret);
  }

  virtual void local_constraints(const ThisType& other, const EntityType& entity,
                                 Constraints::Dirichlet<IntersectionType, RangeFieldType, true>& ret) const
  {
    compute_local_constraints(other, entity, ret);
  }

  virtual void local_constraints(const ThisType& other, const EntityType& entity,
                                 Constraints::Dirichlet<IntersectionType, RangeFieldType, false>& ret) const
  {
    compute_local_constraints(other, entity, ret);
  }

protected:
  mutable Dune::DynamicVector<size_t> tmpMappedRows_;
  mutable Dune::DynamicVector<size_t> tmpMappedCols_;
}; // class ContinuousLagrangeBase


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONTINUOUSLAGRANGE_BASE_HH
