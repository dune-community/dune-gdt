// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
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
#include <dune/xt/grid/boundaryinfo.hh>
#include <dune/xt/grid/intersection.hh>

#include "../interface.hh"

namespace Dune {
namespace GDT {


static constexpr Backends default_cg_backend = default_space_backend;


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

  using typename BaseType::GridLayerType;
  using typename BaseType::EntityType;
  using typename BaseType::PatternType;

private:
  static const constexpr RangeFieldType compare_tolerance_ = 1e-13;

public:
  virtual ~CgSpaceInterface() = default;

  /**
   * \defgroup interface ´´These methods have to be implemented!''
   * @{
   **/
  std::vector<DomainType> lagrange_points(const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp().lagrange_points(entity));
    return this->as_imp().lagrange_points(entity);
  }
  /** @} */

  /**
   * \defgroup provided ´´These methods are provided by the interface for convenience.''
   * @{
   **/
  virtual std::set<size_t> local_dirichlet_DoFs(
      const EntityType& entity,
      const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<GridLayerType>>& boundary_info) const
  {
    std::set<size_t> ret;
    static const XT::Grid::DirichletBoundary dirichlet{};
    const auto lps = lagrange_points(entity);
    const auto intersection_it_end = this->grid_layer().iend(entity);
    for (auto intersection_it = this->grid_layer().ibegin(entity); intersection_it != intersection_it_end;
         ++intersection_it) {
      // only work on dirichlet ones
      const auto& intersection = *intersection_it;
      // actual dirichlet intersections + process boundaries for parallel runs
      if (boundary_info.type(intersection) == dirichlet || (!intersection.neighbor() && !intersection.boundary()))
        for (size_t ii = 0; ii < lps.size(); ++ii) {
          const auto& local_lagrange_point = lps[ii];
          if (XT::Grid::contains(intersection, entity.geometry().global(local_lagrange_point)))
            ret.insert(ii);
        }
    }
    return ret;
  } // ... local_dirichlet_DoFs(...)

  using BaseType::compute_pattern;

  template <class GL, class S, size_t d, size_t r, size_t rC>
  typename std::enable_if<XT::Grid::is_layer<GL>::value, PatternType>::type
  compute_pattern(const GL& grd_layr, const SpaceInterface<S, d, r, rC>& ansatz_space) const
  {
    return BaseType::compute_volume_pattern(grd_layr, ansatz_space);
  }

  using BaseType::local_constraints;

  template <class S, size_t d, size_t r, size_t rC, class ConstraintsType>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*other*/,
                         const EntityType& /*entity*/,
                         ConstraintsType& /*ret*/) const
  {
    static_assert(AlwaysFalse<S>::value, "Not implemented for these constraints!");
  }

  template <class S, size_t d, size_t r, size_t rC>
  void local_constraints(const SpaceInterface<S, d, r, rC>& /*other*/,
                         const EntityType& entity,
                         DirichletConstraints<XT::Grid::extract_intersection_t<GridLayerType>>& ret) const
  {
    const auto local_DoFs = this->local_dirichlet_DoFs(entity, ret.boundary_info());
    if (local_DoFs.size() > 0) {
      const auto global_indices = this->mapper().globalIndices(entity);
      for (const auto& local_DoF : local_DoFs) {
        ret.insert(global_indices[local_DoF]);
      }
    }
  } // ... local_constraints(..., Constraints::Dirichlet<...> ...)
  /** @} */

  static constexpr bool associates_data_with(int codim, std::integral_constant<int, 1>)
  {
    return codim == dimDomain;
  }

  static constexpr bool associates_data_with(int codim, std::integral_constant<int, 2>)
  {
    return dimDomain == 1 ? codim >= 0 : codim > 0;
  }

  static constexpr bool associates_data_with(int codim)
  {
    // actually: return polOrder > 2 ? true : associates_data_with(codim, std::integral_constant<int, polOrder>());
    return codim == 0;
  }

protected:
  std::vector<DomainType> lagrange_points_order_1(const EntityType& entity) const
  {
    // check
    static_assert(polOrder == 1, "Not tested for higher polynomial orders!");
    if (dimRange != 1)
      DUNE_THROW(NotImplemented, "Does not work for higher dimensions");
    assert(this->grid_layer().indexSet().contains(entity));
    // get the basis and reference element
    const auto basis = this->base_function_set(entity);
    typedef typename BaseType::BaseFunctionSetType::RangeType RangeType;
    std::vector<RangeType> tmp_basis_values(basis.size(), RangeType(0));
    const auto& reference_element = ReferenceElements<DomainFieldType, dimDomain>::general(entity.type());
    const auto num_vertices = reference_element.size(dimDomain);
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
}; // class CgSpaceInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CG_INTERFACE_HH
