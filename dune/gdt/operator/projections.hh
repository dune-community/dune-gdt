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
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/assembler/gridwalker.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/space/continuouslagrange.hh>
#include <dune/gdt/space/continuouslagrange/fem.hh>
#include <dune/gdt/space/continuouslagrange/fem-localfunctions.hh>
#include <dune/gdt/space/continuouslagrange/pdelab.hh>
#include <dune/gdt/space/discontinuouslagrange/fem-localfunctions.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace ProjectionOperator {


template <class GridViewImp, class FieldImp = double>
class Lagrange;


template <class GridViewImp, class FieldImp = double>
class LagrangeTraits
{
public:
  typedef Lagrange<GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class LagrangeTraits


/**
 *  \brief  Does a projection using the lagrange points.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridViewImp, class FieldImp>
class Lagrange
{
public:
  typedef LagrangeTraits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

private:
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef FieldVector<DomainFieldType, dimDomain> DomainType;

public:
  Lagrange(const GridViewType& grid_view)
    : grid_view_(grid_view)
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
    // set all dofs to infinity
    const auto infinity = std::numeric_limits<R>::infinity();
    for (size_t ii = 0; ii < range.vector().size(); ++ii)
      range.vector().set_entry(ii, infinity);
    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      const auto lagrange_points = range.space().backend().lagrangePointSet(entity);
      std::vector<DomainType> points(lagrange_points.nop(), DomainType(0));
      for (size_t ii = 0; ii < lagrange_points.nop(); ++ii)
        points[ii] = lagrange_points.point(ii);
      // and do the work (see below)
      apply_local(points, *local_source, local_range_DoF_vector);
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemWrapper< GP, 1, R, r, 1 > ...)

  template <class GP, class R, int r, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, 1>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1>, V>& range) const
  {
    // checks
    typedef ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1> SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // set all dofs to infinity
    const auto infinity = std::numeric_limits<R>::infinity();
    for (size_t ii = 0; ii < range.vector().size(); ++ii)
      range.vector().set_entry(ii, infinity);
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, 1, R, r, 1>, V> RangeFunctionType;
    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity           = *entity_it;
      const auto local_source      = source.local_function(entity);
      auto local_range             = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range.vector();
      const auto lagrange_points   = range.space().lagrange_points(entity);
      // and do the work (see below)
      apply_local(lagrange_points, *local_source, local_range_DoF_vector);
    } // walk the grid
  } // ... apply(... ContinuousLagrangeSpace::FemLocalfunctionsWrapper< GP, 1, R, r, 1 > ...)

private:
  template <class LagrangePointsType, class LocalSourceType, class LocalRangeVectorType>
  void apply_local(const LagrangePointsType& lagrange_points, const LocalSourceType& local_source,
                   LocalRangeVectorType& local_range_DoF_vector) const
  {
    static const unsigned int dimRange = LocalSourceType::dimRange;
    assert(lagrange_points.size() == local_range_DoF_vector.size());
    size_t kk = 0;
    for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
      if (std::isinf(local_range_DoF_vector.get(kk))) {
        const auto& lagrange_point = lagrange_points[ii];
        // evaluate source function
        const auto source_value = local_source.evaluate(lagrange_point);
        for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
          local_range_DoF_vector.set(kk, source_value[jj]);
      } else
        kk += dimRange;
    }
  } // ... apply_local(...)

  const GridViewType& grid_view_;
}; // class Lagrange


template <class GridViewImp, class FieldImp = double>
class L2;


template <class GridViewImp, class FieldImp = double>
class L2Traits
{
public:
  typedef L2<GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class L2Traits


/**
 *  \brief  Does an L2 projection by solving local or global problems.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridViewImp, class FieldImp>
class L2
{
public:
  typedef L2Traits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

private:
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

public:
  L2(const GridViewType& grid_view)
    : grid_view_(grid_view)
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
#if HAVE_EIGEN
    typedef Stuff::LA::EigenDenseMatrix<R> LocalMatrixType;
    typedef Stuff::LA::EigenDenseVector<R> LocalVectorType;
#else // HAVE_EIGEN
    typedef Stuff::LA::CommonDenseMatrix<R> LocalMatrixType;
    typedef Stuff::LA::CommonDenseVector<R> LocalVectorType;
#endif // HAVE_EIGEN
    // clear
    Stuff::Common::clear(range.vector());
    // walk the grid
    RangeType source_value(0);
    std::vector<RangeType> basis_values(range.space().mapper().maxNumDofs());
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      // prepare
      const auto& entity      = *entity_it;
      const auto local_basis  = range.space().base_function_set(entity);
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      LocalMatrixType local_matrix(local_basis.size(), local_basis.size(), R(0));
      LocalVectorType local_vector(local_basis.size(), R(0));
      LocalVectorType local_DoFs(local_basis.size(), R(0));
      // create quadrature
      // guess the polynomial order of the source by hoping that they are the same for all entities
      const size_t integrand_order = std::max(local_source->order(), local_basis.order()) + local_basis.order();
      assert(integrand_order < std::numeric_limits<int>::max());
      const auto& quadrature = QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), int(integrand_order));
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
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix.add_to_entry(
                ii, jj, integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]));
          }
        }
      } // loop over all quadrature points
      // compute local DoFs
      Stuff::LA::Solver<LocalMatrixType>(local_matrix).apply(local_vector, local_DoFs);
      // set local DoFs
      auto local_range_vector = local_range.vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs[ii]);
    } // walk the grid
  } // ... apply(... DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper< GP, p, R, r, 1 > ...)

private:
  const GridViewType& grid_view_;
}; // class L2


template <class GridViewImp, class FieldImp = double>
class Generic;


template <class GridViewImp, class FieldImp = double>
class GenericTraits
{
public:
  typedef Generic<GridViewImp, FieldImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class GenericTraits


/**
 *  \brief  Does a projection by selecting the appropriate Lagrange or L2 operator at compile time.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridViewImp, class FieldImp>
class Generic
{
public:
  typedef GenericTraits<GridViewImp, FieldImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

private:
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

public:
  Generic(const GridViewType& grid_view)
    : lagrange_operator_(grid_view)
    , l2_operator_(grid_view)
  {
  }

  template <class SourceType, class RangeType>
  inline void apply(const SourceType& source, RangeType& range) const
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
  inline void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template <class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V>
  inline void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template <class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V>
  inline void redirect_to_appropriate_operator(
      const Stuff::LocalizableFunctionInterface<E, D, d, RS, rS, rCS>& source,
      DiscreteFunction<DiscontinuousLagrangeSpace::FemLocalfunctionsWrapper<GP, p, RR, rR, rCR>, V>& range) const
  {
    l2_operator_.apply(source, range);
  }

  const Lagrange<GridViewType> lagrange_operator_;
  const L2<GridViewType> l2_operator_;
}; // Generic


// forward, to be used in the traits
template <class GridViewImp, class SourceImp, class RangeImp>
class DirichletLocalizable;


template <class GridViewImp, class SourceImp, class RangeImp>
class DirichletLocalizableTraits
{
  typedef typename RangeImp::SpaceType::Traits T;
  static const unsigned int d = RangeImp::dimDomain;
  typedef typename RangeImp::RangeFieldType R;
  static const unsigned int r  = RangeImp::dimRange;
  static const unsigned int rC = RangeImp::dimRangeCols;
  static_assert(std::is_base_of<ContinuousLagrangeSpaceBase<T, d, R, r, rC>, typename RangeImp::SpaceType>::value,
                "The SpaceType of RangeImp has to be derived from ContinuousLagrangeSpaceBase!");
  static_assert(r == 1, "Not implemeneted for higher dimensions!");
  static_assert(rC == 1, "Not implemeneted for higher dimensions!");
  typedef typename SourceImp::EntityType E;
  typedef typename SourceImp::DomainFieldType D;
  static_assert(SourceImp::dimDomain == d, "Dimensions do not match!");
  static_assert(std::is_same<typename SourceImp::RangeFieldType, R>::value, "Types do not match!");
  static_assert(SourceImp::dimRange == r, "Dimensions do not match!");
  static_assert(SourceImp::dimRangeCols == rC, "Dimensions do not match!");
  static_assert(std::is_base_of<Stuff::LocalizableFunctionInterface<E, D, d, R, r, rC>, SourceImp>::value,
                "SourceImp has to be derived from Stuff::LocalizableFunctionInterface!");

public:
  typedef DirichletLocalizable<GridViewImp, SourceImp, RangeImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef SourceImp SourceType;
  typedef RangeImp RangeType;
}; // class DirichletLocalizableTraits


template <class GridViewImp, class SourceImp, class RangeImp>
class DirichletLocalizable
    : public LocalizableOperatorInterface<DirichletLocalizableTraits<GridViewImp, SourceImp, RangeImp>>,
      public Functor::Codim0<GridViewImp>
{
public:
  typedef DirichletLocalizableTraits<GridViewImp, SourceImp, RangeImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType SourceType;
  typedef typename Traits::RangeType RangeType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;
  typedef Stuff::GridboundaryInterface<IntersectionType> BoundaryInfoType;

public:
  DirichletLocalizable(const GridViewType& grid_view, const BoundaryInfoType& boundary_info, const SourceType& source,
                       RangeType& range)
    : grid_view_(grid_view)
    , boundary_info_(boundary_info)
    , source_(source)
    , range_(range)
  {
  }

  virtual ~DirichletLocalizable()
  {
  }

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE DS_FINAL
  {
    if (entity.hasBoundaryIntersections()) {
      const auto local_dirichlet_DoFs = range_.space().local_dirichlet_DoFs(entity, boundary_info_);
      if (local_dirichlet_DoFs.size() > 0) {
        const auto local_source      = source_.local_function(entity);
        auto local_range             = range_.local_discrete_function(entity);
        auto& local_range_DoF_vector = local_range.vector();
        const auto lagrange_points = range_.space().lagrange_points(entity);
        assert(lagrange_points.size() == local_range_DoF_vector.size());
        for (const size_t& local_DoF_id : local_dirichlet_DoFs)
          local_range_DoF_vector.set(local_DoF_id, local_source->evaluate(lagrange_points[local_DoF_id]));
      }
    }
  } // ... apply_local(...)

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SourceType& source() const
  {
    return source_;
  }

  RangeType& range()
  {
    return range_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  void apply()
  {
    GridWalker<GridViewType> grid_walker(grid_view_);
    grid_walker.add(*this, new ApplyOn::BoundaryEntities<GridViewType>());
    grid_walker.walk();
  }

private:
  const GridViewType& grid_view_;
  const BoundaryInfoType& boundary_info_;
  const SourceType& source_;
  RangeType& range_;
}; // class DirichletLocalizable


template <class GridViewImp>
class Dirichlet;


template <class GridViewImp>
class DirichletTraits
{
public:
  typedef Dirichlet<GridViewImp> derived_type;
  typedef GridViewImp GridViewType;
}; // class DirichletTraits


/**
 *  \brief  Does a dirichlet projection in the sense that the lagrange point set on each entity is matched against
 *          those vertices of the entity which lie on the dirichlet boundary.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template <class GridViewImp>
class Dirichlet
{
public:
  typedef DirichletTraits<GridViewImp> Traits;

  typedef typename Traits::GridViewType GridViewType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;
  typedef Stuff::GridboundaryInterface<typename GridViewType::Intersection> BoundaryInfoType;

public:
  Dirichlet(const GridViewType& grid_view, const BoundaryInfoType& boundary_info)
    : grid_view_(grid_view)
    , boundary_info_(boundary_info)
  {
  }

  template <class R, int r, int rC, class GV, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GV, p, R, r, rC>, V>& range) const
  {
    typedef Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC> SourceType;
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemWrapper<GV, p, R, r, rC>, V> RangeType;
    DirichletLocalizable<GridViewType, SourceType, RangeType> localizable_operator(
        grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  }

  template <class R, int r, int rC, class GV, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& source,
             DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GV, p, R, r, rC>, V>& range) const
  {
    typedef Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC> SourceType;
    typedef DiscreteFunction<ContinuousLagrangeSpace::FemLocalfunctionsWrapper<GV, p, R, r, rC>, V> RangeType;
    DirichletLocalizable<GridViewType, SourceType, RangeType> localizable_operator(
        grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  }

  template <class R, int r, int rC, class GV, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC>& source,
             DiscreteFunction<ContinuousLagrangeSpace::PdelabWrapper<GV, p, R, r, rC>, V>& range) const
  {
    typedef Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, R, r, rC> SourceType;
    typedef DiscreteFunction<ContinuousLagrangeSpace::PdelabWrapper<GV, p, R, r, rC>, V> RangeType;
    DirichletLocalizable<GridViewType, SourceType, RangeType> localizable_operator(
        grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  }

private:
  const GridViewType& grid_view_;
  const BoundaryInfoType& boundary_info_;
}; // class Dirichlet


} // namespace ProjectionOperator
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_PROJECTIONS_HH
