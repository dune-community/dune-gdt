// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_PROJECTIONS_HH
#define DUNE_GDT_OPERATOR_PROJECTIONS_HH

#include <vector>
#include <type_traits>
#include <limits>

#include <dune/common/fvector.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/common/vector.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>

namespace Dune {
namespace GDT {
namespace ProjectionOperator {


/**
 *  \brief  Does a projection using the lagrange points.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridPartImp>
class Lagrange
{
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  Lagrange(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void apply(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
             DiscreteFunction<SpaceInterface<T>, V>& /*range*/) const
  {
    static_assert((Dune::AlwaysFalse<E>::value), "Not implemented for this combination of source and range!");
  }

  template <class GP, class R, int r, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, 1, R, r, 1>, V>& range) const
  {
    // checks
    typedef ContinuousLagrangeSpace::FemWrapper<GP, 1, R, r, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // clear range
    Stuff::Common::clear(range.vector());
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, 1, R, r, 1>, V> RangeFunctionType;
    typedef typename RangeFunctionType::RangeType RangeType;
    RangeType local_source_value(0);
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      const auto lagrange_points = range.space().backend().lagrangePointSet(entity);
      assert(lagrange_points.nop() == local_range_DoF_vector.size());
      // and do the work
      for (size_t ii = 0; ii < lagrange_points.nop(); ++ii) {
        const auto& lagrange_point = lagrange_points.point(ii);
        local_source->evaluate(lagrange_point, local_source_value);
        local_range_DoF_vector.set(ii, local_source_value);
      }
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemWrapper< GP, 1, R, r, 1 > ...)

  template <class GP, class R, int r, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1>, V>& range) const
  {
    // checks
    typedef ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // clear range
    Stuff::Common::clear(range.vector());
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1>, V> RangeFunctionType;
    typedef typename RangeFunctionType::RangeType RangeType;
    RangeType local_source_value(0);
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      const auto lagrange_points = range.space().lagrange_points(entity);
      assert(lagrange_points.size() == local_range_DoF_vector.size());
      // and do the work
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        const auto& lagrange_point = lagrange_points[ii];
        local_source->evaluate(lagrange_point, local_source_value);
        local_range_DoF_vector.set(ii, local_source_value);
      }
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemLocalfunctionsWrapper< GP, 1, R, r, 1 > ...)

private:
  const GridPartType& grid_part_;
}; // class Lagrange


/**
 *  \brief  Does an L2 projection by solving local or global problems.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridPartImp>
class L2
{
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  L2(const GridPartType& grid_part)
    : grid_part_(grid_part)
  {
  }

  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void apply(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
             DiscreteFunction<SpaceInterface<T>, V>& /*range*/) const
  {
    static_assert((Dune::AlwaysFalse<E>::value), "Not implemented for this combination of source and range!");
  }

  template <class GP, int p, class R, int r, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& source,
             DiscreteFunction<DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, R, r, 1>, V>& range) const
  {
    // checks
    typedef DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, R, r, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    typedef typename SpaceType::BaseFunctionSetType::RangeType RangeType;
    // clear
    Stuff::Common::clear(range.vector());
    // walk the grid
    RangeType source_value(0);
    std::vector<RangeType> basis_values(range.space().mapper().maxNumDofs());
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      // prepare
      const auto& entity      = *entity_it;
      const auto local_basis  = range.space().baseFunctionSet(entity);
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      DynamicMatrix<R> local_matrix(local_basis.size(), local_basis.size(), R(0));
      DynamicVector<R> local_vector(local_basis.size(), R(0));
      // create quadrature
      const size_t quadrature_order = std::max(local_source->order(), local_range.order());
      assert((2 * quadrature_order + 1) < std::numeric_limits<int>::max());
      const auto& quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), int(2 * quadrature_order + 1));
      // loop over all quadrature points
      for (const auto& quadrature_point : quadrature) {
        const auto local_point         = quadrature_point.position();
        const auto quadrature_weight   = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        local_source->evaluate(local_point, source_value);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (source_value * basis_values[ii]);
          auto& local_matrix_row = local_matrix[ii];
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix_row[jj] += integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]);
          }
        }
      } // loop over all quadrature points
      // compute local DoFs
      DynamicVector<R> local_DoFs(local_basis.size(), 0);
      local_matrix.solve(local_DoFs, local_vector);
      // set local DoFs
      auto local_range_vector = local_range.vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs[ii]);
    } // walk the grid
  } // ... apply(... DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< GP, p, R, r, 1 > ...)

private:
  const GridPartType& grid_part_;
}; // class L2


/**
 *  \brief  Does a projection by selecting the appropriate Lagrange or L2 operator at compile time.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridPartImp>
class Generic
{
  typedef GridPartImp GridPartType;
  typedef typename GridPartType::template Codim<0>::EntityType EntityType;
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;

public:
  Generic(const GridPartType& grid_part)
    : grid_part_(grid_part)
    , lagrange_operator_(grid_part_)
    , l2_operator_(grid_part_)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range) const
  {
    redirect_to_appropriate_operator(source, range);
  }

private:
  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
                                        DiscreteFunction<SpaceInterface<T>, V>& /*range*/) const
  {
    static_assert((Dune::AlwaysFalse<E>::value),
                  "Could not find an appropriate operator for this combination of source and range!");
  }

  template <class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V>
  void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template <class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V>
  void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template <class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V>
  void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    l2_operator_.apply(source, range);
  }

  const GridPartType& grid_part_;
  const Lagrange<GridPartType> lagrange_operator_;
  const L2<GridPartType> l2_operator_;
}; // Generic


/**
 *  \brief  Does a dirichlet projection in the sense that the lagrange point set on each entity is matched against
 *          those vertices of the entity which lie on the dirichlet boundary.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridPartImp>
class Dirichlet
{
public:
  typedef GridPartImp GridPartType;
  typedef Stuff::GridboundaryInterface<typename GridPartType::IntersectionType> BoundaryInfoType;

private:
  typedef typename GridPartType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridPartType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

  typedef typename GridPartType::template Codim<0>::EntityType EntityType;

public:
  Dirichlet(const GridPartType& grid_part, const BoundaryInfoType& boundary_info)
    : grid_part_(grid_part)
    , boundary_info_(boundary_info)
  {
  }

  template <class E, class D, int d, class R, int r, int rC, class T, class V>
  void apply(const Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>& /*source*/,
             DiscreteFunction<SpaceInterface<T>, V>& /*range*/) const
  {
    static_assert((Dune::AlwaysFalse<E>::value), "Not implemented for this combination of source and range!");
  }

  template <class R, class GP, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, 1, R, 1, 1>, V>& range) const
  {
    // checks
    typedef ContinuousLagrangeSpace::FemWrapper<GP, 1, R, 1, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // clear range
    Stuff::Common::clear(range.vector());
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      // get the lagrange points
      const auto lagrange_points = range.space().backend().lagrangePointSet(entity);
      std::vector<DomainType> points(lagrange_points.nop(), DomainType(0));
      for (size_t ii = 0; ii < lagrange_points.nop(); ++ii)
        points[ii] = lagrange_points.point(ii);
      // and do the work (see below)
      apply_local(entity, points, local_source, local_range_DoF_vector);
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemWrapper< GP, 1, R, 1, 1 > ...) const

  template <class R, class GP, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, 1, 1>, V>& range) const
  {
    // checks
    typedef ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, 1, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // clear range
    Stuff::Common::clear(range.vector());
    // walk the grid
    const auto entity_it_end = grid_part_.template end<0>();
    for (auto entity_it = grid_part_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      const auto lagrange_points   = range.space().lagrange_points(entity);
      // and do the work (see below)
      apply_local(entity, lagrange_points, local_source, local_range_DoF_vector);
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemLocalfunctionsWrapper< GP, 1, R, 1, 1 > ...) const

private:
  template <class LagrangePointsType, class LocalSourceType, class LocalRangeVectorType>
  void apply_local(const EntityType& entity, const LagrangePointsType& lagrange_points,
                   const LocalSourceType& local_source, LocalRangeVectorType& local_range_DoF_vector) const
  {
    assert(lagrange_points.size() == local_range_DoF_vector.size());
    // walk the intersections
    const auto intersection_it_end = grid_part_.iend(entity);
    for (auto intersection_it = grid_part_.ibegin(entity); intersection_it != intersection_it_end; ++intersection_it) {
      const auto& intersection = *intersection_it;
      // only work on boundary intersections
      if (boundary_info_.dirichlet(intersection)) {
        const auto& intersection_geometry = intersection.geometry();
        // and walk its corners (i.e. the vertices in 2d)
        for (int local_intersection_corner_id = 0; local_intersection_corner_id < intersection_geometry.corners();
             ++local_intersection_corner_id) {
          const auto local_vertex = entity.geometry().local(intersection_geometry.corner(local_intersection_corner_id));
          // loop over all local lagrange points
          for (size_t lagrange_point_id = 0; lagrange_point_id < lagrange_points.size(); ++lagrange_point_id) {
            const auto& lagrange_point = lagrange_points[lagrange_point_id];
            // and check for equality
            if (Stuff::Common::FloatCmp::eq(local_vertex, lagrange_point)) {
              // this lagrange point is on the dirichlet boundary
              // so we evaluate the source and set the corresponding local DoF
              local_range_DoF_vector.set(lagrange_point_id, local_source->evaluate(lagrange_point));
            }
          } // loop over all local lagrange points
        } // walk its corners
      } // only work on boundary intersections
    } // walk the intersections
  } // ... apply_local(...)

  const GridPartType& grid_part_;
  const BoundaryInfoType& boundary_info_;
}; // class Dirichlet


} // namespace ProjectionOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROJECTIONS_HH
