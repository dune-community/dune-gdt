// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_HH

#include <vector>
#include <limits>

#include <dune/common/fvector.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/cg/interface.hh>
#include <dune/gdt/spaces/cg/fem.hh>
#include <dune/gdt/spaces/cg/pdelab.hh>
#include <dune/gdt/playground/spaces/dg/fem.hh>
#include <dune/gdt/playground/spaces/dg/pdelab.hh>
#include <dune/gdt/spaces/rt/pdelab.hh>
#include <dune/gdt/spaces/fv/default.hh>
#include <dune/gdt/playground/spaces/block.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


template< class GridViewImp, class FieldImp = double >
class LagrangeProjection;


template< class GridViewImp, class FieldImp = double >
class LagrangeProjectionTraits
{
public:
  typedef LagrangeProjection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class LagrangeProjectionTraits


/**
 *  \brief  Does a projection using the lagrange points.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class LagrangeProjection
{
public:
  typedef LagrangeProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;

private:
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype              DomainFieldType;
  static const unsigned int                         dimDomain = GridViewType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain > DomainType;

public:
  LagrangeProjection(const GridViewType& grid_view)
    : grid_view_(grid_view)
  {}

  template< class E, class D, int d, class R, int r, int rC, class T, class V >
  void apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& /*source*/,
             DiscreteFunction< SpaceInterface< T, d, r, rC >, V >& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse< E >::value, "Not implemented for this combination of source and range!");
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void apply(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::CG::PdelabBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    apply_p(source, range);
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void apply(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::CG::FemBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    apply_p(source, range);
  }

private:
  template< class R, int r, class V, class SpaceType >
  void apply_p(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& source,
             DiscreteFunction< SpaceType, V >& range) const
  {
    // checks
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    // set all dofs to infinity
    const auto infinity = std::numeric_limits< R >::infinity();
    for (size_t ii = 0; ii < range.vector().size(); ++ii)
      range.vector().set_entry(ii, infinity);
    // walk the grid
    const auto entity_it_end = grid_view_.template end< 0 >();
    for (auto entity_it = grid_view_.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      auto& local_range_DoF_vector = local_range->vector();
      const auto lagrange_points = range.space().lagrange_points(entity);
      // and do the work (see below)
      apply_local(lagrange_points, *local_source, local_range_DoF_vector);
    } // walk the grid
  } // ... apply(... Spaces::CG::FemBased< GP, 1, R, r, 1 > ...)

  template< class LagrangePointsType, class LocalSourceType, class LocalRangeVectorType >
  void apply_local(const LagrangePointsType& lagrange_points,
                   const LocalSourceType& local_source,
                   LocalRangeVectorType& local_range_DoF_vector) const
  {
    static const unsigned int dimRange = LocalSourceType::dimRange;
    assert(lagrange_points.size() == local_range_DoF_vector.size());
    size_t kk = 0;
    for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
      if (std::isinf(local_range_DoF_vector[kk])) {
        const auto& lagrange_point = lagrange_points[ii];
        // evaluate source function
        const auto source_value = local_source.evaluate(lagrange_point);
        for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
          local_range_DoF_vector[kk] = source_value[jj];
      }
      else
        kk += dimRange;
    }
  } // ... apply_local(...)

  const GridViewType& grid_view_;
}; // class LagrangeProjection


template< class GridViewImp, class FieldImp = double >
class L2Projection;


template< class GridViewImp, class FieldImp = double >
class L2ProjectionTraits
{
public:
  typedef L2Projection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class L2ProjectionTraits


/**
 *  \brief  Does an L2 projection by solving local or global problems.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class L2Projection
{
public:
  typedef L2ProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;
private:
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridViewType::dimension;

public:
  L2Projection(const GridViewType& grid_view, const size_t over_integrate = 0)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {}

  template< class E, class D, int d, class R, int r, int rC, class T, class V >
  void apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& /*source*/,
             DiscreteFunction< SpaceInterface< T, d, r, rC >, V >& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse< E >::value, "Not implemented for this combination of source and range!");
  }

  template< class GP, int p, class R, int r, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& source,
             DiscreteFunction< Spaces::DG::FemBased< GP, p, R, r, 1 >, V >& range) const
  {
    // checks
    typedef Spaces::DG::FemBased< GP, p, R, r, 1 > SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    apply_local_l2_projection_(source, range);
  } // ... apply(... Spaces::DG::FemBased< ..., 1 > ...)

#if HAVE_DUNE_GRID_MULTISCALE

  template< class GP, int p, class R, int r, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& source,
             DiscreteFunction< Spaces::Block< Spaces::DG::FemBased< GP, p, R, r, 1 > >, V >& range) const
  {
    // checks
    typedef Spaces::Block< Spaces::DG::FemBased< GP, p, R, r, 1 > > SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    apply_local_l2_projection_(source, range);
  }

#endif // HAVE_DUNE_GRID_MULTISCALE

  template< class GP, int p, class R, int r, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, 1 >& source,
             DiscreteFunction< Spaces::DG::PdelabBased< GP, p, R, r, 1 >, V >& range) const
  {
    // checks
    typedef Spaces::DG::PdelabBased< GP, p, R, r, 1 > SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    apply_local_l2_projection_(source, range);
  }

  template< class E, class D, int d, class R, int r, class GV, class V >
  void apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, 1 >& source,
             DiscreteFunction< Spaces::FV::Default< GV, R, r, 1 >, V >& range) const
  {
    typedef Spaces::FV::Default< GV, R, r, 1 > SpaceType;
    static_assert(SpaceType::dimDomain == dimDomain, "Dimensions do not match!");
    apply_local_l2_projection_(source, range);
  } // ... apply(... Spaces::FV::Default< ..., 1 > ...)

  template< class GP, int p, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, FieldType, 1, 1 >& source,
             DiscreteFunction< Spaces::CG::FemBased< GP, p, FieldType, 1, 1 >, V >& range) const
  {
    apply_global_l2_projection_(source, range);
  }

  template< class GP, int p, class V>
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, FieldType, 1, 1 >& source,
             DiscreteFunction< Spaces::CG::PdelabBased< GP, p, FieldType, 1, 1 >, V >& range) const
  {
    apply_global_l2_projection_(source, range);
  }

  template< class GP, int p, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, FieldType, dimDomain, 1 >& source,
             DiscreteFunction< Spaces::RT::PdelabBased< GP, p, FieldType, dimDomain, 1 >, V >& range) const
  {
    apply_global_l2_projection_(source, range);
  } // ... apply(...)

private:
  template< class SourceType, class RangeFunctionType >
  void apply_local_l2_projection_(const SourceType& source, RangeFunctionType& range) const
  {
    typedef typename RangeFunctionType::RangeType RangeType;
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_dense_backend >::MatrixType LocalMatrixType;
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_dense_backend >::VectorType LocalVectorType;
    // clear
    range.vector() *= 0.0;
    // walk the grid
    RangeType source_value(0);
    std::vector< RangeType > basis_values(range.space().mapper().maxNumDofs(), RangeType(0));
    const auto entity_it_end = grid_view_.template end< 0 >();
    for (auto entity_it = grid_view_.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      // prepare
      const auto& entity = *entity_it;
      const auto local_basis = range.space().base_function_set(entity);
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      LocalMatrixType local_matrix(local_basis.size(), local_basis.size(), FieldType(0));
      LocalVectorType local_vector(local_basis.size(), FieldType(0));
      LocalVectorType local_DoFs(local_basis.size(), FieldType(0));
      // create quadrature
      // guess the polynomial order of the source by hoping that they are the same for all entities
      const size_t integrand_order = std::max(local_source->order(), local_basis.order()) + local_basis.order();
      const auto& quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(
            entity.type(), boost::numeric_cast< int >(integrand_order + over_integrate_));
      // loop over all quadrature points
      for (const auto& quadrature_point : quadrature) {
        const auto local_point = quadrature_point.position();
        const auto quadrature_weight = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        local_source->evaluate(local_point, source_value);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (source_value * basis_values[ii]);
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix.add_to_entry(ii,
                                      jj,
                                      integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]));
          }
        }
      } // loop over all quadrature points
      // compute local DoFs
      Stuff::LA::Solver< LocalMatrixType >(local_matrix).apply(local_vector, local_DoFs);
      // set local DoFs
      auto& local_range_vector = local_range->vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector[ii] = local_DoFs[ii];
    } // walk the grid
  } // ... apply_local_l2_projection_(...)

  /**
   * \todo This implementation is not optimal: the MatrixType and VectorType used here do not have to match at all the
   *       VectorType used in RangeFunctionType!
   */
  template< class SourceType, class RangeFunctionType >
  void apply_global_l2_projection_(const SourceType& source, RangeFunctionType& range) const
  {
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_backend >::MatrixType MatrixType;
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_backend >::VectorType VectorType;
    MatrixType lhs(range.space().mapper().size(),
                   range.space().mapper().size(),
                   range.space().compute_volume_pattern());
    VectorType rhs(range.space().mapper().size());

    // walk the grid
    typedef typename RangeFunctionType::LocalfunctionType::RangeType RangeType;
    RangeType source_value(0);
    std::vector< RangeType > basis_values(range.space().mapper().maxNumDofs(), RangeType(0));
    const auto entity_it_end = grid_view_.template end< 0 >();
    for (auto entity_it = grid_view_.template begin< 0 >(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto local_source = source.local_function(entity);
      const auto basis = range.space().base_function_set(entity);
      // do a volume quadrature
      const size_t integrand_order = std::max(std::max(ssize_t(local_source->order()) - 1, ssize_t(0)),
                                              ssize_t(basis.order()))
                                     + basis.order();
      const auto& quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(
            entity.type(), boost::numeric_cast< int >(integrand_order + over_integrate_));
      const auto quadrature_it_end = quadrature.end();
      for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
        const auto xx = quadrature_it->position();
        const FieldType quadrature_weight = quadrature_it->weight();
        const DomainFieldType integration_element = entity.geometry().integrationElement(xx);
        local_source->evaluate(xx, source_value);
        basis.evaluate(xx, basis_values);
        for (size_t ii = 0; ii < basis.size(); ++ii) {
          const size_t global_ii = range.space().mapper().mapToGlobal(entity, ii);
          rhs.add_to_entry(global_ii,
                           integration_element * quadrature_weight * (source_value * basis_values[ii]));
          for (size_t jj = 0; jj < basis.size(); ++jj) {
            const size_t global_jj = range.space().mapper().mapToGlobal(entity, jj);
            lhs.add_to_entry(global_ii,
                             global_jj,
                             integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]));
          }
        }
      } // do a volume quadrature
    } // walk the grid

    // solve
    Stuff::LA::Solver< MatrixType >(lhs).apply(rhs, range.vector());
  } // ... apply_global_l2_projection_(...)

  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class L2Projection


template< class GridViewImp, class FieldImp = double >
class Projection;


template< class GridViewImp, class FieldImp >
class ProjectionTraits
{
public:
  typedef Projection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp FieldType;
}; // class ProjectionTraits


/**
 *  \brief  Does a projection by selecting the appropriate Lagrange or L2 operator at compile time.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class Projection
{
public:
  typedef ProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType    FieldType;
private:
  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype  DomainFieldType;
  static const unsigned int             dimDomain = GridViewType::dimension;

public:
  Projection(const GridViewType& grid_view, const size_t over_integrate = 0)
    : lagrange_operator_(grid_view)
    , l2_operator_(grid_view, over_integrate)
  {}

  template< class SourceType, class RangeType >
  inline void apply(const SourceType& source, RangeType& range) const
  {
    redirect_to_appropriate_operator(source, range);
  }

private:
  template< class E, class D, int d, class R, int r, int rC, class T, class V >
  void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& /*source*/,
                                        DiscreteFunction< SpaceInterface< T, d, r, rC >, V >& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse< E >::value,
                  "Could not find an appropriate operator for this combination of source and range!");
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::CG::FemBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::CG::PdelabBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::DG::FemBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    l2_operator_.apply(source, range);
  }

#if HAVE_DUNE_GRID_MULTISCALE

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::Block< Spaces::DG::FemBased
                                                  < GP, p, RR, rR, rCR > >, V >& range) const
  {
    l2_operator_.apply(source, range);
  }

#endif // HAVE_DUNE_GRID_MULTISCALE

  template< class E, class D, int d, class RS, int rS, int rCS, class GP, int p, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::DG::PdelabBased
                                                  < GP, p, RR, rR, rCR >, V >& range) const
  {
    l2_operator_.apply(source, range);
  }

  template< class E, class D, int d, class RS, int rS, int rCS, class GV, class RR, int rR, int rCR, class V >
  inline void redirect_to_appropriate_operator(const Stuff::LocalizableFunctionInterface< E, D, d, RS, rS, rCS >&
                                                  source,
                                               DiscreteFunction< Spaces::FV::Default< GV, RR, rR, rCR >, V >&
                                                  range) const
  {
    l2_operator_.apply(source, range);
  }

  const LagrangeProjection< GridViewType > lagrange_operator_;
  const L2Projection< GridViewType > l2_operator_;
}; // Projection

template< class SourceType, class RangeType >
void apply_projection(const SourceType& source, RangeType& range) {
  auto& view = range.space().grid_view();
  Projection<typename std::remove_reference<decltype(view)>::type,
      typename RangeType::SpaceType::RangeFieldType>(view).apply(source, range);
}

// forward, to be used in the traits
template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp = double >
class DirichletProjectionLocalizable;


template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class DirichletProjectionLocalizableTraits
{
  typedef typename RangeImp::SpaceType::Traits T;
  static const unsigned int d = RangeImp::dimDomain;
  typedef typename RangeImp::RangeFieldType R;
  static const unsigned int r = RangeImp::dimRange;
  static const unsigned int rC = RangeImp::dimRangeCols;
  static_assert(std::is_base_of< Spaces::CGInterface< T, d, R, r, rC >, typename RangeImp::SpaceType >::value,
                "The SpaceType of RangeImp has to be derived from Spaces::CGInterface!");
  static_assert(r == 1, "Not implemeneted for higher dimensions!");
  static_assert(rC == 1, "Not implemeneted for higher dimensions!");
  typedef typename SourceImp::EntityType E;
  typedef typename SourceImp::DomainFieldType D;
  static_assert(SourceImp::dimDomain == d, "Dimensions do not match!");
  static_assert(std::is_same< typename SourceImp::RangeFieldType, R >::value, "Types do not match!");
  static_assert(SourceImp::dimRange == r, "Dimensions do not match!");
  static_assert(SourceImp::dimRangeCols == rC, "Dimensions do not match!");
  static_assert(std::is_base_of< Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >, SourceImp >::value,
                "SourceImp has to be derived from Stuff::LocalizableFunctionInterface!");
public:
  typedef DirichletProjectionLocalizable< GridViewImp, SourceImp, RangeImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef FieldImp    FieldType;
  typedef SourceImp   SourceType;
  typedef RangeImp    RangeType;
}; // class DirichletProjectionLocalizableTraits


template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class DirichletProjectionLocalizable
  : public LocalizableOperatorInterface< DirichletProjectionLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > >
  , public Stuff::Grid::Functor::Codim0< GridViewImp >
{
public:
  typedef DirichletProjectionLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::RangeType    RangeType;

  typedef typename GridViewType::template Codim< 0 >::Entity      EntityType;
  typedef typename GridViewType::Intersection                     IntersectionType;
  typedef Stuff::Grid::BoundaryInfoInterface< IntersectionType >  BoundaryInfoType;

public:
  DirichletProjectionLocalizable(const GridViewType& grd_vw,
                       const BoundaryInfoType& boundary_info,
                       const SourceType& src,
                       RangeType& rng)
    : grid_view_(grd_vw)
    , boundary_info_(boundary_info)
    , source_(src)
    , range_(rng)
  {}

  virtual ~DirichletProjectionLocalizable() {}

  virtual void apply_local(const EntityType& entity) override final
  {
    if (entity.hasBoundaryIntersections()) {
      const auto local_dirichlet_DoFs = range_.space().local_dirichlet_DoFs(entity, boundary_info_);
      if (local_dirichlet_DoFs.size() > 0) {
        const auto local_source = source_.local_function(entity);
        auto local_range = range_.local_discrete_function(entity);
        auto& local_range_DoF_vector = local_range->vector();
        const auto lagrange_points = range_.space().lagrange_points(entity);
        assert(lagrange_points.size() == local_range_DoF_vector.size());
        for (const size_t& local_DoF_id : local_dirichlet_DoFs)
          local_range_DoF_vector[local_DoF_id] = local_source->evaluate(lagrange_points[local_DoF_id]);
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
    Stuff::Grid::Walker< GridViewType > grid_walker(grid_view_);
    grid_walker.add(*this, new Stuff::Grid::ApplyOn::BoundaryEntities< GridViewType >());
    grid_walker.walk();
  }

private:
  const GridViewType& grid_view_;
  const BoundaryInfoType& boundary_info_;
  const SourceType& source_;
  RangeType& range_;
}; // class DirichletProjectionLocalizable


template< class GridViewImp >
class DirichletProjection;


template< class GridViewImp >
class DirichletProjectionTraits
{
public:
  typedef DirichletProjection< GridViewImp > derived_type;
  typedef GridViewImp GridViewType;
}; // class DirichletProjectionTraits


/**
 *  \brief  Does a dirichlet projection in the sense that the lagrange point set on each entity is matched against
 *          those vertices of the entity which lie on the dirichlet boundary.
 *  \note   This use of the lagrange points is known to fail for polynomial orders higher than 1.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp >
class DirichletProjection
{
public:
  typedef DirichletProjectionTraits< GridViewImp > Traits;

  typedef typename Traits::GridViewType GridViewType;

  typedef typename GridViewType::template Codim< 0 >::Entity                  EntityType;
  typedef typename GridViewType::ctype                                        DomainFieldType;
  static const unsigned int                                                   dimDomain = GridViewType::dimension;
  typedef Stuff::Grid::BoundaryInfoInterface< typename GridViewType::Intersection > BoundaryInfoType;

public:
  DirichletProjection(const GridViewType& grid_view, const BoundaryInfoType& boundary_info)
    : grid_view_(grid_view)
    , boundary_info_(boundary_info)
  {}

  template< class R, int r, int rC, class GV, int p, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
             DiscreteFunction< Spaces::CG::FemBased< GV, p, R, r, rC >, V >& range) const
  {
    typedef Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC > SourceType;
    typedef DiscreteFunction< Spaces::CG::FemBased< GV, p, R, r, rC >, V >           RangeType;
    DirichletProjectionLocalizable< GridViewType, SourceType, RangeType >
        localizable_operator(grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  }

  template< class R, int r, int rC, class GV, int p, class V >
  void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
             DiscreteFunction< Spaces::CG::PdelabBased< GV, p, R, r, rC >, V >& range) const
  {
    typedef Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC > SourceType;
    typedef DiscreteFunction< Spaces::CG::PdelabBased< GV, p, R, r, rC >, V >        RangeType;
    DirichletProjectionLocalizable< GridViewType, SourceType, RangeType >
        localizable_operator(grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  }

private:
  const GridViewType& grid_view_;
  const BoundaryInfoType& boundary_info_;
}; // class DirichletProjection

template< class SourceType, class RangeSpaceType, class V >
void apply_dirichlet_projection(
    const DSG::BoundaryInfoInterface< typename RangeSpaceType::GridViewType::Intersection>& boundary_info,
    const SourceType& source, DiscreteFunction<RangeSpaceType, V>& range) {
  auto& view = range.space().grid_view();
  DirichletProjection<typename std::remove_reference<decltype(*view)>::type>(*view, boundary_info).apply(source, range);
}

} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_HH
