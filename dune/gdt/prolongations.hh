// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_PROLONGATIONS_HH
#define DUNE_GDT_OPERATORS_PROLONGATIONS_HH

#include <vector>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>

#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/grid/intersection.hh>
#include <dune/stuff/la/container.hh>
#include <dune/stuff/la/solver.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/grid/search.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/cg/fem.hh>
#include <dune/gdt/spaces/cg/pdelab.hh>
#include <dune/gdt/spaces/fv/defaultproduct.hh>


namespace Dune {
namespace GDT {

namespace Spaces {
namespace DG {
template <class GridPartImp, int polynomialOrder, class RangeFieldImp, size_t rangeDim, size_t rangeDimCols>
class FemBased;
}
template <class SpaceImp>
class Block;
}

namespace Operators {


/**
 *  \note The automatic detection of the right integration order might fail, so you might want to specify
 *        over_integrate. The reason is that in order to locally evaluate the source we first have to create a
 *        quadrature, the correct order of wich we guess by taking the sources order on the first entity.
 *  \note We would have liked to do something like this and match on implementations of SpaceInterface:\code
template< class T, class VS, class GPR, int pR, class RR, size_t rR, size_t rCR, class VR >
void apply(const ConstDiscreteFunction< SpaceInterface< T >, VS >& source,
           DiscreteFunction< Spaces::DG::FemBased< GPR, pR, RR, rR, rCR >, VR >& range) const
{
  static_assert(Dune::AlwaysFalse< T >::value, "Not implemented for this combination of source and range!");
}\endcode
 *        but that gave compile errors (the compiler just could not match the first argument for whatever reason). This
 *        is why we need all combinations of spaces below which are just compile time checks and forwards.
 *
 *  \todo refactor like projections, drop note above
 */
template <class GridViewType>
class L2Prolongation
{
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

public:
  L2Prolongation(const GridViewType& grid_view)
    : grid_view_(grid_view)
  {
  }

  // Source: Spaces::CG::FemBased
  // Range:  Spaces::DG::FemBased

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  void apply(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, RS, rS, rCS>, VS>& /*source*/,
             DiscreteFunction<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>, VR>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<GPS>::value, "Not implemented for this combination of source and range!");
  }

  template <class GPS, int pS, class R, size_t r, size_t rC, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, R, r, rC>, VS>& source,
                    DiscreteFunction<Spaces::DG::FemBased<GPR, pR, R, r, rC>, VR>& range) const
  {
    prolong_onto_dg_fem_localfunctions_wrapper(source, range);
  }

  // Source: Spaces::DG::FemBased
  // Range:  Spaces::DG::FemBased

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  void apply(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>, VS>& /*source*/,
             DiscreteFunction<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>, VR>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<GPS>::value, "Not implemented for this combination of source and range!");
  }

  template <class GPS, int pS, class R, size_t r, size_t rC, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, R, r, rC>, VS>& source,
                    DiscreteFunction<Spaces::DG::FemBased<GPR, pR, R, r, rC>, VR>& range) const
  {
    prolong_onto_dg_fem_localfunctions_wrapper(source, range);
  }

  template <class GPS, int pS, class R, size_t r, size_t rC, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPS, pS, R, r, rC>>, VS>& source,
                    DiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPR, pR, R, r, rC>>, VR>& range) const
  {
    prolong_onto_dg_fem_localfunctions_wrapper(source, range);
  }

  template <class GPS, int pS, class R, size_t r, size_t rC, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPS, pS, R, r, rC>>, VS>& source,
                    DiscreteFunction<Spaces::DG::FemBased<GPR, pR, R, r, rC>, VR>& range) const
  {
    prolong_onto_dg_fem_localfunctions_wrapper(source, range);
  }

  template <class GVS, class RS, size_t rS, size_t rCS, class VS, class GVR, class RR, size_t rR, size_t rCR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::FV::DefaultProduct<GVS, RS, rS, rCS>, VS>& source,
                    DiscreteFunction<Spaces::FV::DefaultProduct<GVR, RR, rR, rCR>, VR>& range) const
  {
    prolong_onto_fv(source, range);
  }

private:
  template <class SourceFunctionType, class RangeFunctionType>
  void prolong_onto_dg_fem_localfunctions_wrapper(const SourceFunctionType& source, RangeFunctionType& range) const
  {
    typedef typename RangeFunctionType::DomainType DomainType;
    typedef typename RangeFunctionType::RangeType RangeType;
    typedef typename RangeFunctionType::RangeFieldType RangeFieldType;
    typedef typename Stuff::LA::Container<RangeFieldType, Stuff::LA::default_dense_backend>::MatrixType LocalMatrixType;
    typedef typename Stuff::LA::Container<RangeFieldType, Stuff::LA::default_dense_backend>::VectorType LocalVectorType;
    // clear
    range.vector() *= 0.0;
    // create search in the source grid part
    typedef typename SourceFunctionType::SpaceType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch<SourceGridViewType> EntitySearch;
    EntitySearch entity_search(source.space().grid_view());
    // guess the polynomial order of the source by hoping that they are the same for all entities
    const size_t source_order = source.local_function(*source.space().grid_view().template begin<0>())->order();
    // walk the grid
    RangeType source_value(0);
    std::vector<RangeType> basis_values(range.space().mapper().maxNumDofs());
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      // prepare
      const auto& entity     = *entity_it;
      const auto local_basis = range.space().base_function_set(entity);
      auto local_range = range.local_discrete_function(entity);
      LocalMatrixType local_matrix(local_basis.size(), local_basis.size(), RangeFieldType(0));
      LocalVectorType local_vector(local_basis.size(), RangeFieldType(0));
      LocalVectorType local_DoFs(local_basis.size(), RangeFieldType(0));
      // create quadrature
      const auto integrand_order = std::max(source_order, local_basis.order()) + local_basis.order();
      const auto& quadrature =
          QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), boost::numeric_cast<int>(integrand_order));
      // get global quadrature points
      std::vector<DomainType> quadrature_points;
      for (const auto& quadrature_point : quadrature)
        quadrature_points.emplace_back(entity.geometry().global(quadrature_point.position()));
      // get source entities
      const auto source_entity_ptr_unique_ptrs = entity_search(quadrature_points);
      assert(source_entity_ptr_unique_ptrs.size() >= quadrature_points.size());
      // loop over all quadrature points
      size_t pp = 0;
      for (const auto& quadrature_point : quadrature) {
        const auto local_point         = quadrature_point.position();
        const auto quadrature_weight   = quadrature_point.weight();
        const auto integration_element = entity.geometry().integrationElement(local_point);
        // evaluate source
        const auto& source_entity_ptr_unique_ptr = source_entity_ptr_unique_ptrs[pp];
        if (source_entity_ptr_unique_ptr) {
          const auto source_entity_ptr = *source_entity_ptr_unique_ptr;
          const auto& source_entity    = *source_entity_ptr;
          const auto local_source = source.local_function(source_entity);
          local_source->evaluate(source_entity.geometry().local(entity.geometry().global(local_point)), source_value);
        } else
          source_value *= 0.0;
        // evaluate
        local_basis.evaluate(local_point, basis_values);
        // compute integrals
        for (size_t ii = 0; ii < local_basis.size(); ++ii) {
          local_vector[ii] += integration_element * quadrature_weight * (source_value * basis_values[ii]);
          for (size_t jj = 0; jj < local_basis.size(); ++jj) {
            local_matrix.add_to_entry(
                ii, jj, integration_element * quadrature_weight * (basis_values[ii] * basis_values[jj]));
          }
        }
        ++pp;
      } // loop over all quadrature points
      // compute local DoFs
      try {
        Stuff::LA::Solver<LocalMatrixType>(local_matrix).apply(local_vector, local_DoFs);
      } catch (Stuff::Exceptions::linear_solver_failed& ee) {
        DUNE_THROW(Exceptions::prolongation_error,
                   "L2 prolongation failed because a local matrix could not be inverted!\n\n"
                       << "This was the original error: "
                       << ee.what());
      }
      // set local DoFs
      auto local_range_vector = local_range->vector();
      assert(local_range_vector.size() == local_DoFs.size());
      for (size_t ii = 0; ii < local_range_vector.size(); ++ii)
        local_range_vector.set(ii, local_DoFs.get_entry(ii));
    } // walk the grid
  } // ... prolong_onto_dg_fem_localfunctions_wrapper(...)

  template <class SourceFunctionType, class RangeFunctionType>
  void prolong_onto_fv(const SourceFunctionType& source, RangeFunctionType& range) const
  {
    typedef typename RangeFunctionType::DomainType DomainType;
    typedef typename RangeFunctionType::RangeType RangeType;
    // create search in the source grid part
    typedef typename SourceFunctionType::SpaceType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch<SourceGridViewType> EntitySearch;
    EntitySearch entity_search(source.space().grid_view());
    // walk the grid
    RangeType source_value(0);
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      auto local_range   = range.local_discrete_function(entity);
      // get global quadrature points
      std::vector<DomainType> quadrature_points(1, entity.geometry().center());
      // get source entities
      const auto source_entity_ptr_unique_ptrs = entity_search(quadrature_points);
      assert(source_entity_ptr_unique_ptrs.size() >= 1);
      const auto& source_entity_unique_ptr = source_entity_ptr_unique_ptrs[0];
      if (source_entity_unique_ptr) {
        const auto source_entity_ptr = *source_entity_unique_ptr;
        const auto local_source = source.local_function(*source_entity_ptr);
        local_source->evaluate(source_entity_ptr->geometry().local(entity.geometry().center()), source_value);
      } else
        source_value *= 0.0;
      // set local DoFs
      auto local_range_vector = local_range->vector();
      assert(local_range_vector.size() == RangeFunctionType::dimRange);
      for (size_t ii = 0; ii < RangeFunctionType::dimRange; ++ii)
        local_range_vector.set(ii, source_value[ii]);
    } // walk the grid
  } // ... prolong_onto_fv(...)

  const GridViewType& grid_view_;
}; // class L2Prolongation


/**
 *  \note We would have liked to do something like this and match on implementations of SpaceInterface:\code
template< class T, class VS, class GPR, int pR, class RR, size_t rR, size_t rCR, class VR >
void apply(const ConstDiscreteFunction< SpaceInterface< T >, VS >& source,
           DiscreteFunction< Spaces::CG::FemBased< GPR, pR, RR, rR, rCR >, VR >& range) const
{
  static_assert(Dune::AlwaysFalse< T >::value, "Not implemented for this combination of source and range!");
}\endcode
 *        but that gave compile errors (the compiler just could not match the first argument for whatever reason). This
 *        is why we need all combinations of spaces below which are just compile time checks and forwards.
 *
 *  \todo refactor like projections, drop note above
 */
template <class GridViewType>
class LagrangeProlongation
{
public:
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

  LagrangeProlongation(const GridViewType& grid_view)
    : grid_view_(grid_view)
  {
  }

  // Source: Spaces::CG::FemBased
  // Range:  Spaces::CG::FemBased

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  void apply(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, RS, rS, rCS>, VS>& /*source*/,
             DiscreteFunction<Spaces::CG::FemBased<GPR, pR, RR, rR, rCR>, VR>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<GPS>::value, "Not implemented for this combination of source and range!");
  }

  template <class GPS, int pS, class R, size_t r, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, R, r, 1>, VS>& source,
                    DiscreteFunction<Spaces::CG::FemBased<GPR, pR, R, r, 1>, VR>& range) const
  {
    redirect_to_appropriate_apply(source, range);
  }

  // Source: Spaces::DG::FemBased
  // Range:  Spaces::CG::FemBased

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  void apply(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>, VS>& /*source*/,
             DiscreteFunction<Spaces::CG::FemBased<GPR, pR, RR, rR, rCR>, VR>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<GPS>::value, "Not implemented for this combination of source and range!");
  }

  template <class GPS, int pS, class R, size_t r, class VS, class GPR, int pR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, R, r, 1>, VS>& source,
                    DiscreteFunction<Spaces::CG::FemBased<GPR, pR, R, r, 1>, VR>& range) const
  {
    redirect_to_appropriate_apply(source, range);
  }

  // Source: Spaces::CG::PdelabBased
  // Range:  Spaces::CG::PdelabBased

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  void apply(const ConstDiscreteFunction<Spaces::CG::PdelabBased<GPS, pS, RS, rS, rCS>, VS>& /*source*/,
             DiscreteFunction<Spaces::CG::PdelabBased<GPR, pR, RR, rR, rCR>, VR>& /*range*/) const
  {
    static_assert(Dune::AlwaysFalse<GPS>::value, "Not implemented for this combination of source and range!");
  }

  template <class GPS, int pS, class R, size_t r, size_t rC, class VS, class GPR, class VR>
  inline void apply(const ConstDiscreteFunction<Spaces::CG::PdelabBased<GPS, pS, R, r, rC>, VS>& source,
                    DiscreteFunction<Spaces::CG::PdelabBased<GPR, 1, R, r, rC>, VR>& range) const
  {
    redirect_to_appropriate_apply(source, range);
  }

private:
  template <class SourceType, class RangeType>
  void redirect_to_appropriate_apply(const SourceType& source, RangeType& range) const
  {
    // create search in the source grid part
    typedef typename SourceType::SpaceType::GridViewType SourceGridViewType;
    typedef Stuff::Grid::EntityInlevelSearch<SourceGridViewType> EntitySearch;
    EntitySearch entity_search(source.space().grid_view());
    // set all range dofs to infinity
    const auto infinity = std::numeric_limits<typename RangeType::RangeFieldType>::infinity();
    for (size_t ii = 0; ii < range.vector().size(); ++ii)
      range.vector().set_entry(ii, infinity);
    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get global lagrange point coordinates
      const auto lagrange_point_set = range.space().lagrange_points(entity);
      typedef FieldVector<typename SourceGridViewType::ctype, SourceGridViewType::dimension> DomainType;
      std::vector<DomainType> lagrange_points(lagrange_point_set.size());
      for (size_t ii = 0; ii < lagrange_point_set.size(); ++ii)
        lagrange_points[ii] = entity.geometry().global(lagrange_point_set[ii]);
      // get source entities
      const auto source_entity_ptrs = entity_search(lagrange_points);
      assert(source_entity_ptrs.size() == lagrange_points.size());
      // get range
      auto local_range            = range.local_discrete_function(entity);
      auto local_range_DoF_vector = local_range->vector();
      // do the actual work (see below)
      apply_local(source, lagrange_points, source_entity_ptrs, local_range_DoF_vector);
    } // walk the grid
  } // ... redirect_to_appropriate_apply(...)

  template <class SourceType, class LagrangePointsType, class EntityPointers, class LocalDoFVectorType>
  void apply_local(const SourceType& source, const LagrangePointsType& lagrange_points,
                   const EntityPointers& source_entity_ptr_unique_ptrs, LocalDoFVectorType& range_DoF_vector) const
  {
    static const size_t dimRange = SourceType::dimRange;
    size_t kk = 0;
    assert(source_entity_ptr_unique_ptrs.size() >= lagrange_points.size());
    for (size_t ii = 0; ii < lagrange_points.size(); ++ii) {
      if (std::isinf(range_DoF_vector.get(kk))) {
        const auto& global_point = lagrange_points[ii];
        // evaluate source function
        const auto& source_entity_ptr_unique_ptr = source_entity_ptr_unique_ptrs[ii];
        if (source_entity_ptr_unique_ptr) {
          const auto source_entity_ptr  = *source_entity_ptr_unique_ptr;
          const auto& source_entity     = *source_entity_ptr;
          const auto local_source_point = source_entity.geometry().local(global_point);
          const auto local_source       = source.local_function(source_entity);
          const auto source_value = local_source->evaluate(local_source_point);
          for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
            range_DoF_vector.set(kk, source_value[jj]);
        } else
          for (size_t jj = 0; jj < dimRange; ++jj, ++kk)
            range_DoF_vector.set(kk, 0.0);
      } else
        kk += dimRange;
    }
  } // ... apply_local(...)

  const GridViewType& grid_view_;
}; // class LagrangeProlongation


template <class GridViewType>
class Prolongation
{
public:
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

  Prolongation(const GridViewType& grid_view)
    : l2_prolongation_operator_(grid_view)
    , lagrange_prolongation_operator_(grid_view)
  {
  }

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range) const
  {
    redirect_to_appropriate_operator(source, range);
  }

private:
  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void
  redirect_to_appropriate_operator(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, RS, rS, rCS>, VS>& source,
                                   DiscreteFunction<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    l2_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void
  redirect_to_appropriate_operator(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>, VS>& source,
                                   DiscreteFunction<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    l2_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void redirect_to_appropriate_operator(
      const ConstDiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>>, VS>& source,
      DiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>>, VR>& range) const
  {
    l2_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void redirect_to_appropriate_operator(
      const ConstDiscreteFunction<Spaces::Block<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>>, VS>& source,
      DiscreteFunction<Spaces::DG::FemBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    l2_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void
  redirect_to_appropriate_operator(const ConstDiscreteFunction<Spaces::CG::FemBased<GPS, pS, RS, rS, rCS>, VS>& source,
                                   DiscreteFunction<Spaces::CG::FemBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    lagrange_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void
  redirect_to_appropriate_operator(const ConstDiscreteFunction<Spaces::DG::FemBased<GPS, pS, RS, rS, rCS>, VS>& source,
                                   DiscreteFunction<Spaces::CG::FemBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    lagrange_prolongation_operator_.apply(source, range);
  }

  template <class GPS, int pS, class RS, size_t rS, size_t rCS, class VS, class GPR, int pR, class RR, size_t rR,
            size_t rCR, class VR>
  inline void redirect_to_appropriate_operator(
      const ConstDiscreteFunction<Spaces::CG::PdelabBased<GPS, pS, RS, rS, rCS>, VS>& source,
      DiscreteFunction<Spaces::CG::PdelabBased<GPR, pR, RR, rR, rCR>, VR>& range) const
  {
    lagrange_prolongation_operator_.apply(source, range);
  }

  template <class GVS, class RS, size_t rS, size_t rCS, class VS, class GVR, class RR, size_t rR, size_t rCR, class VR>
  inline void redirect_to_appropriate_operator(
      const ConstDiscreteFunction<Spaces::FV::DefaultProduct<GVS, RS, rS, rCS>, VS>& source,
      DiscreteFunction<Spaces::FV::DefaultProduct<GVR, RR, rR, rCR>, VR>& range) const
  {
    l2_prolongation_operator_.apply(source, range);
  }

  const L2Prolongation<GridViewType> l2_prolongation_operator_;
  const LagrangeProlongation<GridViewType> lagrange_prolongation_operator_;
}; // class Prolongation


template <class GridViewType, class SourceType, class RangeType>
void prolong(const GridViewType& grid_view, const SourceType& source, RangeType& range)
{
  const Prolongation<GridViewType> prolongation_operator(grid_view);
  prolongation_operator.apply(source, range);
}


template <class SourceType, class RangeType>
void prolong(const SourceType& source, RangeType& range)
{
  const Prolongation<typename RangeType::SpaceType::GridViewType> prolongation_operator(range.space().grid_view());
  prolongation_operator.apply(source, range);
}


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_PROLONGATIONS_HH
