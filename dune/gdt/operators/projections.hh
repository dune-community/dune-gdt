// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROJECTIONS_HH
#define DUNE_GDT_OPERATORS_PROJECTIONS_HH

#include <vector>
#include <limits>

#if HAVE_TBB
# include <tbb/blocked_range.h>
# include <tbb/parallel_reduce.h>
# include <tbb/tbb_stddef.h>
#endif

#include <dune/common/fvector.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/cg/interface.hh>
#include <dune/gdt/spaces/dg/interface.hh>
#include <dune/gdt/spaces/fv/interface.hh>
#include <dune/gdt/spaces/rt/interface.hh>
#include <dune/gdt/playground/spaces/block.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forwards
template< class GridViewImp, class FieldImp = double >
class LagrangeProjection;

template< class GridViewImp, class FieldImp = double >
class L2Projection;

template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp = double >
class DirichletProjectionLocalizable;

template< class GridViewImp >
class DirichletProjection;

template< class GridViewImp, class FieldImp = double >
class Projection;


namespace internal {


template< class GridViewImp, class FieldImp = double >
class LagrangeProjectionTraits
{
public:
  typedef LagrangeProjection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp                                 GridViewType;
  typedef FieldImp                                    FieldType;
};


template< class GridViewImp, class FieldImp = double >
class L2ProjectionTraits
{
public:
  typedef L2Projection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp                           GridViewType;
  typedef FieldImp                              FieldType;
};


template< class GridViewImp, class FieldImp >
class ProjectionTraits
{
public:
  typedef Projection< GridViewImp, FieldImp > derived_type;
  typedef GridViewImp                         GridViewType;
  typedef FieldImp                            FieldType;
};


template< class GridViewImp >
class DirichletProjectionTraits
{
public:
  typedef DirichletProjection< GridViewImp > derived_type;
  typedef GridViewImp                        GridViewType;
};



template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class DirichletProjectionLocalizableTraits
{
  typedef typename RangeImp::SpaceType::Traits T;
  static const size_t d = RangeImp::dimDomain;
  typedef typename RangeImp::RangeFieldType R;
  static const size_t r = RangeImp::dimRange;
  static const size_t rC = RangeImp::dimRangeCols;
  static_assert(is_cg_space< typename RangeImp::SpaceType >::value,
                "The SpaceType of RangeImp has to be derived from Spaces::CGInterface!");
  static_assert(r == 1, "Not implemeneted for higher dimensions!");
  static_assert(rC == 1, "Not implemeneted for higher dimensions!");
  typedef typename SourceImp::EntityType E;
  typedef typename SourceImp::DomainFieldType D;
  static_assert(SourceImp::dimDomain == d, "Dimensions do not match!");
  static_assert(std::is_same< typename SourceImp::RangeFieldType, R >::value, "Types do not match!");
  static_assert(SourceImp::dimRange == r, "Dimensions do not match!");
  static_assert(SourceImp::dimRangeCols == rC, "Dimensions do not match!");
  static_assert(Stuff::is_localizable_function< SourceImp >::value,
                "SourceImp has to be derived from Stuff::LocalizableFunctionInterface!");
public:
  typedef DirichletProjectionLocalizable< GridViewImp, SourceImp, RangeImp > derived_type;
  typedef GridViewImp                                                        GridViewType;
  typedef FieldImp                                                           FieldType;
  typedef SourceImp                                                          SourceType;
  typedef RangeImp                                                           RangeType;
}; // class DirichletProjectionLocalizableTraits


} // namespace internal


/**
 *  \brief  Does a projection using the lagrange points.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class LagrangeProjection
{
public:
  typedef internal::LagrangeProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType                               GridViewType;
  typedef typename Traits::FieldType                                  FieldType;
private:
  typedef typename GridViewType::template Codim< 0 >::Entity          EntityType;
  typedef typename GridViewType::ctype                                DomainFieldType;
  static const size_t                                                 dimDomain = GridViewType::dimension;
  typedef FieldVector< DomainFieldType, dimDomain >                   DomainType;

public:
  LagrangeProjection(const GridViewType& grid_view)
    : grid_view_(grid_view)
  {}

  template< class R, size_t r, size_t rC, class S, class V >
  inline void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
                    DiscreteFunction< S, V >& range) const
  {
    redirect_apply(range.space(), source, range);
  }

private:
  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::CGInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                      const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
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
      assert(lagrange_points.size() == local_range_DoF_vector.size());
      size_t kk = 0;
      for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
        if (std::isinf(local_range_DoF_vector.get(kk))) {
          const auto& lagrange_point = lagrange_points[ii];
          // evaluate source function
          const auto source_value = local_source->evaluate(lagrange_point);
          for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
            local_range_DoF_vector.set(kk, source_value[jj]);
        } else
          kk += dimRange;
      }
    } // walk the grid
  } // ... redirect_apply(...)

  const GridViewType& grid_view_;
}; // class LagrangeProjection


/**
 *  \brief  Does an L2 projection by solving local or global problems.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class L2Projection
{
public:
  typedef L2Projection< GridViewImp, FieldImp >        ThisType;
  typedef internal::L2ProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType                         GridViewType;
  typedef typename Traits::FieldType                            FieldType;
private:
  typedef typename GridViewType::template Codim< 0 >::Entity    EntityType;
  typedef typename GridViewType::ctype                          DomainFieldType;
  static const size_t                                           dimDomain = GridViewType::dimension;

public:
  L2Projection(const GridViewType& grid_view, const size_t over_integrate = 0)
    : grid_view_(grid_view)
    , over_integrate_(over_integrate)
  {}

  /**
   * \brief Applies an L2 projection.
   *
   *        Extracts the space type and calls the correct method, \sa redirect_apply.
   */
  template< class R, size_t r, size_t rC, class S, class V >
  inline void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
                    DiscreteFunction< S, V >& range) const
  {
    redirect_apply(range.space(), source, range);
  }

  template< class R, size_t r, size_t rC, class S, class V >
  inline void apply(const Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain, EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
                    DiscreteFunction< S, V >& range) const
  {
    redirect_apply(range.space(), source, range);
  }

private:
  template< class T, class R, size_t dimRange, class S, class V >
  inline void redirect_apply(const Spaces::CGInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                             const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                             DiscreteFunction< S, V >& range) const
  {
    apply_global_l2_projection(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::DGInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                      const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::FVInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                      const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::FVInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                      const Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain, EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection_expression_checkerboard(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::ProductFVInterface< T >& /*space*/,
                      const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::ProductFVInterface< T >& /*space*/,
                      const Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain, EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection_expression_checkerboard_fv(source, range);
  }

  template< class T, class R, size_t dimRange, class S, class V >
  void redirect_apply(const Spaces::RTInterface< T, dimDomain, dimRange, 1 >& /*space*/,
                      const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                      DiscreteFunction< S, V >& range) const
  {
    apply_local_l2_projection(source, range);
  }

#if HAVE_DUNE_GRID_MULTISCALE

  /**
   * \brief Extracts the local space type and redirects to the correct method.
   */
  template< class L, class R, size_t dimRange, class S, class V >
  inline void redirect_apply(const Spaces::Block< L >& space,
                             const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& source,
                             DiscreteFunction< S, V >& range) const
  {
    redirect_apply(*(space.local_spaces().at(0)), source, range);
  }

#endif // HAVE_DUNE_GRID_MULTISCALE

private:
#if HAVE_TBB
  template< class SourceType, class RangeFunctionType, class PartitioningType >
  struct Body
  {
    Body(const ThisType& projection_operator,
         const PartitioningType& partitioning,
         const SourceType& source,
         RangeFunctionType& range)
      : projection_operator_(projection_operator)
      , partitioning_(partitioning)
      , source_(source)
      , range_(range)
    {}

    Body(Body& other, tbb::split /*split*/)
      : projection_operator_(other.projection_operator_)
      , partitioning_(other.partitioning_)
      , source_(other.source_)
      , range_(other.range_)
    {}

    void operator()(const tbb::blocked_range< std::size_t > &range) const
    {
      // for all partitions in tbb-range
      for(std::size_t p = range.begin(); p != range.end(); ++p) {
        auto partition = partitioning_.partition(p);
        projection_operator_.walk_grid_parallel(source_, range_, partition);
      }
    }

    void join(Body& /*other*/)
    {}

    const ThisType& projection_operator_;
    const PartitioningType& partitioning_;
    const SourceType& source_;
    RangeFunctionType& range_;
  }; // struct Body

  template< class SourceType, class RangeFunctionType, class PartitioningType >
  struct BodyFV
  {
    BodyFV(const ThisType& projection_operator,
         const PartitioningType& partitioning,
         const SourceType& source,
         RangeFunctionType& range)
      : projection_operator_(projection_operator)
      , partitioning_(partitioning)
      , source_(source)
      , range_(range)
    {}

    BodyFV(BodyFV& other, tbb::split /*split*/)
      : projection_operator_(other.projection_operator_)
      , partitioning_(other.partitioning_)
      , source_(other.source_)
      , range_(other.range_)
    {}

    void operator()(const tbb::blocked_range< std::size_t > &range) const
    {
      // for all partitions in tbb-range
      for(std::size_t p = range.begin(); p != range.end(); ++p) {
        auto partition = partitioning_.partition(p);
        projection_operator_.walk_grid_parallel_expression_checkerboard_fv(source_, range_, partition);
      }
    }

    void join(BodyFV& /*other*/)
    {}

    const ThisType& projection_operator_;
    const PartitioningType& partitioning_;
    const SourceType& source_;
    RangeFunctionType& range_;
  }; // struct BodyFV
#endif //HAVE_TBB

  template< class SourceType, class RangeFunctionType, class EntityRange >
  void walk_grid_parallel(const SourceType& source, RangeFunctionType& range, const EntityRange& entity_range) const
  {
    typedef typename RangeFunctionType::RangeType RangeType;
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_dense_backend >::MatrixType LocalMatrixType;
    typedef typename Stuff::LA::Container< FieldType, Stuff::LA::default_dense_backend >::VectorType LocalVectorType;
    RangeType source_value(0);
    std::vector< RangeType > basis_values(range.space().mapper().maxNumDofs(), RangeType(0));
#ifdef __INTEL_COMPILER
    const auto it_end = entity_range.end();
    for (auto it = entity_range.begin(); it != it_end; ++it) {
      const EntityType& entity = *it;
#else
    for (const EntityType& entity : entity_range) {
#endif
      // prepare
      const auto local_basis = range.space().base_function_set(entity);
      const auto local_source = source.local_function(entity);
      auto local_range = range.local_discrete_function(entity);
      LocalMatrixType local_matrix(local_basis.size(), local_basis.size(), FieldType(0));
      LocalVectorType local_vector(local_basis.size(), FieldType(0));
      LocalVectorType local_DoFs(local_basis.size(), FieldType(0));
      // create quadrature
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
      try {
        Stuff::LA::Solver< LocalMatrixType >(local_matrix).apply(local_vector, local_DoFs);
      } catch (Stuff::Exceptions::linear_solver_failed& ee) {
        DUNE_THROW(Exceptions::projection_error,
                   "L2 projection failed because a local matrix could not be inverted!\n\n"
                   << "This was the original error: " << ee.what());
      }
      // set local DoFs
      auto local_range_vector = local_range->vector();
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs[ii]);
    } // walk the grid
  } // void walk_grid_parallel

  template< class SourceType, class RangeFunctionType, class EntityRange >
  void walk_grid_parallel_expression_checkerboard_fv(const SourceType& source, RangeFunctionType& range, const EntityRange& entity_range) const
  {
    typedef typename RangeFunctionType::RangeType RangeType;
    std::shared_ptr< const typename SourceType::RangeType > source_expression;
    RangeType source_value(0);
#ifdef __INTEL_COMPILER
    const auto it_end = entity_range.end();
    for (auto it = entity_range.begin(); it != it_end; ++it) {
      const EntityType& entity = *it;
#else
    for (const EntityType& entity : entity_range) {
#endif
      // prepare
      const auto local_source = source.local_global_function(entity);
      auto local_range = range.local_discrete_function(entity);
      RangeType local_vector(0);
      // create quadrature
      const size_t integrand_order = local_source->order();
      const auto& quadrature = QuadratureRules< DomainFieldType, dimDomain >::rule(
            entity.type(), boost::numeric_cast< int >(integrand_order));
      // loop over all quadrature points
      for (const auto& quadrature_point : quadrature) {
        const auto local_point = quadrature_point.position();
        const auto quadrature_weight = quadrature_point.weight();
      // evaluate
      local_source->evaluate(local_point, source_expression);
      // compute integrals
      source_expression->evaluate(entity.geometry().global(local_point), source_value);
      local_vector.axpy(quadrature_weight, source_value);
      }
      // set local DoFs
      auto local_range_vector = local_range->vector();
      assert(local_range_vector.size() == RangeFunctionType::dimRange);
      for (size_t ii = 0; ii < RangeFunctionType::dimRange; ++ii)
        local_range_vector.set(ii, local_vector[ii]);
    } // walk the grid
  } // void walk_grid_parallel_expression_checkerboard_fv

  template< class SourceType, class RangeFunctionType >
  void apply_local_l2_projection(const SourceType& source, RangeFunctionType& range) const
  {
    // clear
    std::fill(range.vector().begin(), range.vector().end(), 0.0);
#if HAVE_TBB
    // create partitioning
    const auto num_partitions = DSC_CONFIG_GET("threading.partition_factor", 1u)
                                * DS::threadManager().current_threads();
    RangedPartitioning< GridViewType, 0 > partitioning(range.space().grid_view(), num_partitions);
    tbb::blocked_range< std::size_t > blocked_range(0, partitioning.partitions());
    Body< SourceType, RangeFunctionType, RangedPartitioning< GridViewType, 0 > > body(*this, partitioning, source, range);
   // walk the grid
    tbb::parallel_reduce(blocked_range, body);
#else // HAVE_TBB
    walk_grid_parallel(source, range, GridViewType::elements(range.space().grid_view()));
#endif // HAVE_TBB
  } // ... apply_local_l2_projection(...)

  template< class SourceType, class RangeFunctionType >
  void apply_local_l2_projection_expression_checkerboard_fv(const SourceType& source, RangeFunctionType& range) const
  {
    // clear
    std::fill(range.vector().begin(), range.vector().end(), 0.0);
#if HAVE_TBB
    // create partitioning
    const auto num_partitions = DSC_CONFIG_GET("threading.partition_factor", 1u)
                                * DS::threadManager().current_threads();
    RangedPartitioning< GridViewType, 0 > partitioning(range.space().grid_view(), num_partitions);
    tbb::blocked_range< std::size_t > blocked_range(0, partitioning.partitions());
    BodyFV< SourceType, RangeFunctionType, RangedPartitioning< GridViewType, 0 > > body(*this, partitioning, source, range);
    // walk the grid
    tbb::parallel_reduce(blocked_range, body);
#else // HAVE_TBB
    walk_grid_parallel_expression_checkerboard_fv(source, range, GridViewType::elements(range.space().grid_view()));
#endif // HAVE_TBB
  } // ... apply_local_l2_projection_expression_checkerboard(...)

  template< class SourceType, class RangeFunctionType >
  void apply_global_l2_projection(const SourceType& source, RangeFunctionType& range) const
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
    try {
      Stuff::LA::Solver< MatrixType >(lhs).apply(rhs, range.vector());
    } catch (Stuff::Exceptions::linear_solver_failed& ee) {
      DUNE_THROW(Exceptions::projection_error,
                 "L2 projection failed because the projection matrix could not be inverted!\n\n"
                 << "This was the original error: " << ee.what());
    }
  } // ... apply_global_l2_projection(...)

  const GridViewType& grid_view_;
  const size_t over_integrate_;
}; // class L2Projection


/**
 *  \brief  Does a projection by selecting the appropriate Lagrange or L2 operator at compile time.
 *  \note   If you add other dimension/polorder/space combinations, do not forget to add a testcase in
 *          tests/operators.cc!
 */
template< class GridViewImp, class FieldImp >
class Projection
{
public:
  typedef internal::ProjectionTraits< GridViewImp, FieldImp > Traits;
  typedef typename Traits::GridViewType                       GridViewType;
  typedef typename Traits::FieldType                          FieldType;
private:
  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const size_t                                         dimDomain = GridViewType::dimension;

public:
  Projection(const GridViewType& grid_view, const size_t over_integrate = 0)
    : lagrange_operator_(grid_view)
    , l2_operator_(grid_view, over_integrate)
  {}

  template< class R, size_t r, size_t rC, class S, class V >
  inline void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
                    DiscreteFunction< S, V >& range) const
  {
    redirect_apply(range.space(), source, range);
  }

  template< class R, size_t r, size_t rC, class S, class V >
  inline void apply(const Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain, EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
                    DiscreteFunction< S, V >& range) const
  {
    redirect_apply(range.space(), source, range);
  }

private:
  template< class T, class R, size_t dimRange, size_t dimRangeCols, class S, class V >
  inline void redirect_apply(const Spaces::CGInterface< T, dimDomain, dimRange, dimRangeCols >& /*space*/,
                             const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, dimRangeCols >& source,
                             DiscreteFunction< S, V >& range) const
  {
    lagrange_operator_.apply(source, range);
  }

  template< class T, class R, size_t dimRange, size_t dimRangeCols, class S, class V >
  inline void redirect_apply(const SpaceInterface< T, dimDomain, dimRange, dimRangeCols >& /*space*/,
                             const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, dimRange, dimRangeCols >& source,
                             DiscreteFunction< S, V >& range) const
  {
    l2_operator_.apply(source, range);
  }

  template< class T, class R, size_t dimRange, size_t dimRangeCols, class S, class V >
  inline void redirect_apply(const SpaceInterface< T, dimDomain, dimRange, dimRangeCols >& /*space*/,
                             const Stuff::GlobalFunctionValuedFunctionInterface< EntityType, DomainFieldType, dimDomain, EntityType, DomainFieldType, dimDomain, R, dimRange, dimRangeCols >& source,
                             DiscreteFunction< S, V >& range) const
  {
    l2_operator_.apply(source, range);
  }

  const LagrangeProjection< GridViewType > lagrange_operator_;
  const L2Projection< GridViewType > l2_operator_;
}; // Projection


template< class GridViewImp, class SourceImp, class RangeImp, class FieldImp >
class DirichletProjectionLocalizable
  : public LocalizableOperatorInterface
        < internal::DirichletProjectionLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > >
  , public Stuff::Grid::Functor::Codim0< GridViewImp >
{
public:
  typedef internal::DirichletProjectionLocalizableTraits< GridViewImp, SourceImp, RangeImp, FieldImp > Traits;
  typedef typename Traits::GridViewType                                                                GridViewType;
  typedef typename Traits::SourceType                                                                  SourceType;
  typedef typename Traits::RangeType                                                                   RangeType;
  typedef typename Stuff::Grid::Functor::Codim0< GridViewImp >::EntityType                             EntityType;
  typedef Stuff::Grid::BoundaryInfoInterface< typename GridViewType::Intersection >                    BoundaryInfoType;

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
  typedef internal::DirichletProjectionTraits< GridViewImp >                        Traits;
  typedef typename Traits::GridViewType                                             GridViewType;
  typedef typename GridViewType::template Codim< 0 >::Entity                        EntityType;
  typedef typename GridViewType::ctype                                              DomainFieldType;
  static const size_t                                                               dimDomain = GridViewType::dimension;
  typedef Stuff::Grid::BoundaryInfoInterface< typename GridViewType::Intersection > BoundaryInfoType;

public:
  DirichletProjection(const GridViewType& grid_view, const BoundaryInfoType& boundary_info)
    : grid_view_(grid_view)
    , boundary_info_(boundary_info)
  {}

  template< class R, size_t r, size_t rC, class S, class V >
   void apply(const Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC >& source,
              DiscreteFunction< S, V >& range) const
  {
    typedef Stuff::LocalizableFunctionInterface< EntityType, DomainFieldType, dimDomain, R, r, rC > SourceType;
    typedef DiscreteFunction< S, V >                                                                RangeType;
    DirichletProjectionLocalizable< GridViewType, SourceType, RangeType >
        localizable_operator(grid_view_, boundary_info_, source, range);
    localizable_operator.apply();
  } // ... apply(...)

private:
  const GridViewType& grid_view_;
  const BoundaryInfoType& boundary_info_;
}; // class DirichletProjection


} // namespace Operators


template< class E, class D, size_t d, class R, size_t r, size_t rC, class S, class V >
inline void project_lagrange(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
                             DiscreteFunction< S, V >& range)
{
  Operators::LagrangeProjection< typename S::GridViewType, R >(range.space().grid_view()).apply(source, range);
}


template< class E, class D, size_t d, class R, size_t r, size_t rC, class S, class V >
inline void project_l2(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
                       DiscreteFunction< S, V >& range)
{
  Operators::L2Projection< typename S::GridViewType, R >(range.space().grid_view()).apply(source, range);
}


template< class E, class D, size_t d, class R, size_t r, size_t rC, class S, class V >
inline void project(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
                    DiscreteFunction< S, V >& range)
{
  Operators::Projection< typename S::GridViewType, R >(range.space().grid_view()).apply(source, range);
}

template< class E, class D, size_t d, class RE, class RD, size_t Rd, class R, size_t r, size_t rC, class S, class V >
inline void project(const Stuff::GlobalFunctionValuedFunctionInterface< E, D, d, RE, RD, Rd, R, r, rC >& source,
                    DiscreteFunction< S, V >& range)
{
  Operators::Projection< typename S::GridViewType, R >(range.space().grid_view()).apply(source, range);
}


template< class E, class D, size_t d, class R, size_t r, size_t rC, class S, class V >
inline void project_dirichlet(const DSG::BoundaryInfoInterface< typename S::GridViewType::Intersection >& boundary_info,
                              const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
                              DiscreteFunction< S, V >& range)
{
  Operators::DirichletProjection< typename S::GridViewType >(range.space().grid_view(), boundary_info).apply(source, range);
}


namespace Operators {


template< class SourceType, class RangeType >
void
  DUNE_DEPRECATED_MSG("Use Dune::GDT::project() instead (08.02.2015)! (drop the Operators NS!)")
     apply_projection(const SourceType& source, RangeType& range)
{
  auto& view = range.space().grid_view();
  Projection< typename std::remove_reference<decltype(view)>::type,
      typename RangeType::SpaceType::RangeFieldType>(view).apply(source, range);
}


template< class SourceType, class RangeSpaceType, class V >
void
  DUNE_DEPRECATED_MSG("Use Dune::GDT::project_dirichlet() instead (08.02.2015)! (drop the Operators NS!)")
     apply_dirichlet_projection(
    const DSG::BoundaryInfoInterface< typename RangeSpaceType::GridViewType::Intersection >& boundary_info,
    const SourceType& source,
    DiscreteFunction<RangeSpaceType, V >& range)
{
  auto& view = range.space().grid_view();
  DirichletProjection< typename std::remove_reference< decltype(view) >::type >(view, boundary_info).apply(source, range);
}


template< class GV, class S, class R >
    DirichletProjectionLocalizable< GV, S, R >
make_localizable_dirichlet_projection(const GV& grid_view,
                                      const Stuff::Grid::BoundaryInfoInterface< typename GV::Intersection >& boundary_info,
                                      const S& source,
                                      R& range)
{
  return DirichletProjectionLocalizable< GV, S, R >(grid_view, boundary_info, source, range);
}


template< class GV >
    DirichletProjection< GV >
make_dirichlet_projection(const GV& grid_view,
                          const Stuff::Grid::BoundaryInfoInterface< typename GV::Intersection >& boundary_info)
{
  return DirichletProjection< GV >(grid_view, boundary_info);
}


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROJECTIONS_HH
