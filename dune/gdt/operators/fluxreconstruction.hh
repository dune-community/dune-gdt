// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_OPERATORS_RECONSTRUCTIONS_HH
#define DUNE_GDT_OPERATORS_RECONSTRUCTIONS_HH

#include <type_traits>
#include <limits>

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/functions/constant.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/integrands/elliptic-ipdg.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

namespace Dune {
namespace GDT {


template <class GridLayerType, class DiffusionFactorType, class DiffusionTensorType = void>
class DiffusiveFluxReconstructionOperator;


template <class GridLayerType, class DiffusionFactorType, class DiffusionTensorType>
class DiffusiveFluxReconstructionOperator
{
  static_assert(GridLayerType::dimension == 2, "Only implemented for dimDomain 2 at the moment!");
  static_assert(Dune::XT::Functions::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be tagged as XT::Functions::is_localizable_function!");
  static_assert(Dune::XT::Functions::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be tagged as XT::Functions::is_localizable_function!");

public:
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename DiffusionFactorType::RangeFieldType FieldType;
  typedef typename DiffusionFactorType::DomainType DomainType;

private:
  static_assert(dimDomain == 2, "Not implemented!");

public:
  DiffusiveFluxReconstructionOperator(const GridLayerType& grid_layer,
                                      const DiffusionFactorType& diffusion_factor,
                                      const DiffusionTensorType& diffusion_tensor,
                                      const size_t over_integrate = 0)
    : grid_layer_(grid_layer)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , over_integrate_(over_integrate)
  {
  }

  template <class GL, class V>
  void
  apply(const XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1>& source,
        DiscreteFunction<DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>, V>& range) const
  {
    const auto& rtn0_space = range.space();
    auto& range_vector = range.vector();
    const FieldType infinity = std::numeric_limits<FieldType>::infinity();
    for (size_t ii = 0; ii < range_vector.size(); ++ii)
      range_vector[ii] = infinity;
    const LocalEllipticIpdgIntegrands::Inner<DiffusionFactorType, DiffusionTensorType> inner_evaluation(
        diffusion_factor_, diffusion_tensor_);
    const LocalEllipticIpdgIntegrands::BoundaryLHS<DiffusionFactorType, DiffusionTensorType> boundary_evaluation(
        diffusion_factor_, diffusion_tensor_);
    const XT::Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, FieldType, 1> constant_one(1);
    DomainType normal(0);
    DomainType xx_entity(0);
    DynamicMatrix<FieldType> tmp_matrix(1, 1, 0);
    DynamicMatrix<FieldType> tmp_matrix_en_en(1, 1, 0);
    DynamicMatrix<FieldType> tmp_matrix_en_ne(1, 1, 0);
    std::vector<typename DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>::BaseFunctionSetType::RangeType>
        basis_values(rtn0_space.mapper().maxNumDofs(),
                     typename DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>::BaseFunctionSetType::RangeType(0));
    // walk the grid
    const auto entity_it_end = grid_layer_.template end<0>();
    for (auto entity_it = grid_layer_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto local_DoF_indices = rtn0_space.local_DoF_indices(entity);
      const auto global_DoF_indices = rtn0_space.mapper().globalIndices(entity);
      assert(global_DoF_indices.size() == local_DoF_indices.size());
      const auto local_diffusion_factor = diffusion_factor_.local_function(entity);
      const auto local_diffusion_tensor = diffusion_tensor_.local_function(entity);
      const auto local_source = source.local_function(entity);
      const auto local_basis = rtn0_space.base_function_set(entity);
      const auto local_constant_one = constant_one.local_function(entity);
      // walk the intersections
      const auto intersection_it_end = grid_layer_.iend(entity);
      for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbor = intersection.outside();
          if (grid_layer_.indexSet().index(entity) < grid_layer_.indexSet().index(neighbor)) {
            const auto local_diffusion_factor_neighbor = diffusion_factor_.local_function(neighbor);
            const auto local_diffusion_tensor_neighbor = diffusion_tensor_.local_function(neighbor);
            const auto local_source_neighbor = source.local_function(neighbor);
            const auto local_constant_one_neighbor = constant_one.local_function(neighbor);
            const size_t local_intersection_index = intersection.indexInInside();
            const size_t local_DoF_index = local_DoF_indices[local_intersection_index];
            // do a face quadrature
            FieldType lhs = 0;
            FieldType rhs = 0;
            const size_t integrand_order = inner_evaluation.order(*local_diffusion_factor,
                                                                  *local_diffusion_tensor,
                                                                  *local_diffusion_factor_neighbor,
                                                                  *local_diffusion_tensor_neighbor,
                                                                  *local_constant_one,
                                                                  *local_source,
                                                                  *local_constant_one_neighbor,
                                                                  *local_source_neighbor);
            const auto& quadrature = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
                intersection.type(), boost::numeric_cast<int>(integrand_order + over_integrate_));
            const auto quadrature_it_end = quadrature.end();
            for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
              const auto& xx_intersection = quadrature_it->position();
              xx_entity = intersection.geometryInInside().global(xx_intersection);
              normal = intersection.unitOuterNormal(xx_intersection);
              const FieldType integration_factor = intersection.geometry().integrationElement(xx_intersection);
              const FieldType weigth = quadrature_it->weight();
              // evalaute
              local_basis.evaluate(xx_entity, basis_values);
              const auto& basis_value = basis_values[local_DoF_index];
              tmp_matrix_en_en *= 0.0;
              tmp_matrix_en_ne *= 0.0;
              inner_evaluation.evaluate(*local_diffusion_factor,
                                        *local_diffusion_tensor,
                                        *local_diffusion_factor_neighbor,
                                        *local_diffusion_tensor_neighbor,
                                        *local_constant_one,
                                        *local_source,
                                        *local_constant_one_neighbor,
                                        *local_source_neighbor,
                                        intersection,
                                        xx_intersection,
                                        tmp_matrix_en_en, // <- we are interested in this one
                                        tmp_matrix,
                                        tmp_matrix_en_ne, // <- and this one
                                        tmp_matrix);
              // compute integrals
              assert(tmp_matrix_en_en.rows() >= 1);
              assert(tmp_matrix_en_en.cols() >= 1);
              assert(tmp_matrix_en_ne.rows() >= 1);
              assert(tmp_matrix_en_ne.cols() >= 1);
              lhs += integration_factor * weigth * (basis_value * normal);
              rhs += integration_factor * weigth * (tmp_matrix_en_en[0][0] + tmp_matrix_en_ne[0][0]);
            } // do a face quadrature
            // set DoF
            const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
            // and make sure we are the first to do so
            assert(!(range_vector[global_DoF_index] < infinity));
            range_vector[global_DoF_index] = rhs / lhs;
          }
        } else if (intersection.boundary() && !intersection.neighbor()) {
          const size_t local_intersection_index = intersection.indexInInside();
          const size_t local_DoF_index = local_DoF_indices[local_intersection_index];
          // do a face quadrature
          FieldType lhs = 0;
          FieldType rhs = 0;
          const size_t integrand_order = boundary_evaluation.order(
              *local_diffusion_factor, *local_diffusion_tensor, *local_source, *local_constant_one);
          const auto& quadrature = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
              intersection.type(), boost::numeric_cast<int>(integrand_order + over_integrate_));
          const auto quadrature_it_end = quadrature.end();
          for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
            const auto xx_intersection = quadrature_it->position();
            normal = intersection.unitOuterNormal(xx_intersection);
            const FieldType integration_factor = intersection.geometry().integrationElement(xx_intersection);
            const FieldType weigth = quadrature_it->weight();
            xx_entity = intersection.geometryInInside().global(xx_intersection);
            // evalaute
            local_basis.evaluate(xx_entity, basis_values);
            const auto& basis_value = basis_values[local_DoF_index];
            tmp_matrix *= 0.0;
            boundary_evaluation.evaluate(*local_diffusion_factor,
                                         *local_diffusion_tensor,
                                         *local_constant_one,
                                         *local_source,
                                         intersection,
                                         xx_intersection,
                                         tmp_matrix);
            // compute integrals
            assert(tmp_matrix.rows() >= 1);
            assert(tmp_matrix.cols() >= 1);
            lhs += integration_factor * weigth * (basis_value * normal);
            rhs += integration_factor * weigth * tmp_matrix[0][0];
          } // do a face quadrature
          // set DoF
          const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
          // and make sure we are the first to do so
          assert(!(range_vector[global_DoF_index] < infinity));
          range_vector[global_DoF_index] = rhs / lhs;
        } else
          DUNE_THROW(XT::Common::Exceptions::internal_error, "Unknown intersection type!");
      } // walk the intersections
    } // walk the grid
  } // ... apply(...)

private:
  const GridLayerType& grid_layer_;
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const size_t over_integrate_;
}; // class DiffusiveFluxReconstructionOperator


/**
 *  \todo Add more static checks that GridLayerType and LocalizableFunctionType match.
 *  \todo Derive from operator interfaces.
 */
template <class GridLayerType, class LocalizableFunctionType>
class DiffusiveFluxReconstructionOperator<GridLayerType, LocalizableFunctionType, void>
{
  static_assert(GridLayerType::dimension == 2, "Only implemented for dimDomain 2 at the moment!");
  static_assert(XT::Functions::is_localizable_function<LocalizableFunctionType>::value,
                "LocalizableFunctionType has to be tagged as XT::Functions::is_localizable_function!");

public:
  using EntityType = XT::Grid::extract_entity_t<GridLayerType>;
  typedef typename GridLayerType::ctype DomainFieldType;
  static const size_t dimDomain = GridLayerType::dimension;
  typedef typename LocalizableFunctionType::RangeFieldType FieldType;
  typedef typename LocalizableFunctionType::DomainType DomainType;

private:
  static_assert(dimDomain == 2, "Not implemented!");

public:
  DiffusiveFluxReconstructionOperator(const GridLayerType& grid_layer,
                                      const LocalizableFunctionType& diffusion,
                                      const size_t over_integrate = 0)
    : grid_layer_(grid_layer)
    , diffusion_(diffusion)
    , over_integrate_(over_integrate)
  {
  }

  template <class GL, class V>
  void
  apply(const XT::Functions::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1>& source,
        DiscreteFunction<DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>, V>& range) const
  {
    const auto& rtn0_space = range.space();
    auto& range_vector = range.vector();
    const FieldType infinity = std::numeric_limits<FieldType>::infinity();
    for (size_t ii = 0; ii < range_vector.size(); ++ii)
      range_vector[ii] = infinity;
    const LocalEllipticIpdgIntegrands::Inner<LocalizableFunctionType> inner_evaluation(diffusion_);
    const LocalEllipticIpdgIntegrands::BoundaryLHS<LocalizableFunctionType> boundary_evaluation(diffusion_);
    const XT::Functions::ConstantFunction<EntityType, DomainFieldType, dimDomain, FieldType, 1> constant_one(1);
    DomainType normal(0);
    DomainType xx_entity(0);
    DynamicMatrix<FieldType> tmp_matrix(1, 1, 0);
    DynamicMatrix<FieldType> tmp_matrix_en_en(1, 1, 0);
    DynamicMatrix<FieldType> tmp_matrix_en_ne(1, 1, 0);
    std::vector<typename DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>::BaseFunctionSetType::RangeType>
        basis_values(rtn0_space.mapper().maxNumDofs(),
                     typename DunePdelabRtSpaceWrapper<GL, 0, FieldType, dimDomain>::BaseFunctionSetType::RangeType(0));
    // walk the grid
    const auto entity_it_end = grid_layer_.template end<0>();
    for (auto entity_it = grid_layer_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      const auto local_DoF_indices = rtn0_space.local_DoF_indices(entity);
      const auto global_DoF_indices = rtn0_space.mapper().globalIndices(entity);
      assert(global_DoF_indices.size() == local_DoF_indices.size());
      const auto local_diffusion = diffusion_.local_function(entity);
      const auto local_source = source.local_function(entity);
      const auto local_basis = rtn0_space.base_function_set(entity);
      const auto local_constant_one = constant_one.local_function(entity);
      // walk the intersections
      const auto intersection_it_end = grid_layer_.iend(entity);
      for (auto intersection_it = grid_layer_.ibegin(entity); intersection_it != intersection_it_end;
           ++intersection_it) {
        const auto& intersection = *intersection_it;
        if (intersection.neighbor() && !intersection.boundary()) {
          const auto neighbor = intersection.outside();
          if (grid_layer_.indexSet().index(entity) < grid_layer_.indexSet().index(neighbor)) {
            const auto local_diffusion_neighbor = diffusion_.local_function(neighbor);
            const auto local_source_neighbor = source.local_function(neighbor);
            const auto local_constant_one_neighbor = constant_one.local_function(neighbor);
            const size_t local_intersection_index = intersection.indexInInside();
            const size_t local_DoF_index = local_DoF_indices[local_intersection_index];
            // do a face quadrature
            FieldType lhs = 0;
            FieldType rhs = 0;
            const size_t integrand_order = inner_evaluation.order(*local_diffusion,
                                                                  *local_diffusion_neighbor,
                                                                  *local_constant_one,
                                                                  *local_source,
                                                                  *local_constant_one_neighbor,
                                                                  *local_source_neighbor);
            const auto& quadrature = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
                intersection.type(), boost::numeric_cast<int>(integrand_order + over_integrate_));
            const auto quadrature_it_end = quadrature.end();
            for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
              const auto& xx_intersection = quadrature_it->position();
              xx_entity = intersection.geometryInInside().global(xx_intersection);
              normal = intersection.unitOuterNormal(xx_intersection);
              const auto integration_factor = intersection.geometry().integrationElement(xx_intersection);
              const auto weigth = quadrature_it->weight();
              // evalaute
              local_basis.evaluate(xx_entity, basis_values);
              const auto& basis_value = basis_values[local_DoF_index];
              tmp_matrix_en_en *= 0.0;
              tmp_matrix_en_ne *= 0.0;
              inner_evaluation.evaluate(*local_diffusion,
                                        *local_diffusion_neighbor,
                                        *local_constant_one,
                                        *local_source,
                                        *local_constant_one_neighbor,
                                        *local_source_neighbor,
                                        intersection,
                                        xx_intersection,
                                        tmp_matrix_en_en, // <- we are interested in this one
                                        tmp_matrix,
                                        tmp_matrix_en_ne, // <- and this one
                                        tmp_matrix);
              // compute integrals
              assert(tmp_matrix_en_en.rows() >= 1);
              assert(tmp_matrix_en_en.cols() >= 1);
              assert(tmp_matrix_en_ne.rows() >= 1);
              assert(tmp_matrix_en_ne.cols() >= 1);
              lhs += integration_factor * weigth * (basis_value * normal);
              rhs += integration_factor * weigth * (tmp_matrix_en_en[0][0] + tmp_matrix_en_ne[0][0]);
            } // do a face quadrature
            // set DoF
            const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
            // and make sure we are the first to do so
            assert(!(range_vector[global_DoF_index] < infinity));
            range_vector[global_DoF_index] = rhs / lhs;
          }
        } else if (intersection.boundary() && !intersection.neighbor()) {
          const size_t local_intersection_index = intersection.indexInInside();
          const size_t local_DoF_index = local_DoF_indices[local_intersection_index];
          // do a face quadrature
          FieldType lhs = 0;
          FieldType rhs = 0;
          const size_t integrand_order =
              boundary_evaluation.order(*local_diffusion, *local_source, *local_constant_one);
          const auto& quadrature = QuadratureRules<DomainFieldType, dimDomain - 1>::rule(
              intersection.type(), boost::numeric_cast<int>(integrand_order + over_integrate_));
          const auto quadrature_it_end = quadrature.end();
          for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
            const auto xx_intersection = quadrature_it->position();
            normal = intersection.unitOuterNormal(xx_intersection);
            const auto integration_factor = intersection.geometry().integrationElement(xx_intersection);
            const auto weigth = quadrature_it->weight();
            xx_entity = intersection.geometryInInside().global(xx_intersection);
            // evalaute
            local_basis.evaluate(xx_entity, basis_values);
            const auto& basis_value = basis_values[local_DoF_index];
            tmp_matrix *= 0.0;
            boundary_evaluation.evaluate(
                *local_diffusion, *local_constant_one, *local_source, intersection, xx_intersection, tmp_matrix);
            // compute integrals
            assert(tmp_matrix.rows() >= 1);
            assert(tmp_matrix.cols() >= 1);
            lhs += integration_factor * weigth * (basis_value * normal);
            rhs += integration_factor * weigth * tmp_matrix[0][0];
          } // do a face quadrature
          // set DoF
          const size_t global_DoF_index = global_DoF_indices[local_DoF_index];
          // and make sure we are the first to do so
          assert(!(range_vector[global_DoF_index] < infinity));
          range_vector[global_DoF_index] = rhs / lhs;
        } else
          DUNE_THROW(XT::Common::Exceptions::internal_error, "Unknown intersection type!");
      } // walk the intersections
    } // walk the grid
  } // ... apply(...)

private:
  const GridLayerType& grid_layer_;
  const LocalizableFunctionType& diffusion_;
  const size_t over_integrate_;
}; // class DiffusiveFluxReconstructionOperator


template <class GL, class DF, class DT>
DiffusiveFluxReconstructionOperator<GL, DF, DT> make_diffusive_flux_reconstruction_operator(
    const GL& grid_layer, const DF& diffusion_factor, const DT& diffusion_tensor, const size_t over_integrate = 0)
{
  return DiffusiveFluxReconstructionOperator<GL, DF, DT>(
      grid_layer, diffusion_factor, diffusion_tensor, over_integrate);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_RECONSTRUCTIONS_HH
