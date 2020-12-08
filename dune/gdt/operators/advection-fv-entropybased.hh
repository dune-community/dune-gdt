// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
#define DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH

#include <boost/config.hpp>

#include <dune/xt/common/mkl.hh>
#include <dune/xt/common/numeric.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


template <class OperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesOperator
  : public OperatorInterface<typename OperatorImp::AGV,
                             OperatorImp::s_r,
                             OperatorImp::s_rC,
                             OperatorImp::r_r,
                             OperatorImp::r_rC,
                             typename OperatorImp::F,
                             typename OperatorImp::M,
                             typename OperatorImp::SGV,
                             typename OperatorImp::RGV>
{
public:
  using ThisType = EntropicCoordinatesOperator;
  using BaseType = OperatorInterface<typename OperatorImp::AGV,
                                     OperatorImp::s_r,
                                     OperatorImp::s_rC,
                                     OperatorImp::r_r,
                                     OperatorImp::r_rC,
                                     typename OperatorImp::F,
                                     typename OperatorImp::M,
                                     typename OperatorImp::SGV,
                                     typename OperatorImp::RGV>;

  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using OperatorType = OperatorImp;
  using InverseHessianOperatorType = InverseHessianOperatorImp;

  EntropicCoordinatesOperator(const OperatorType& opertr, const InverseHessianOperatorType& inverse_hessian_operator)
    : BaseType(opertr.parameter_type() + inverse_hessian_operator.parameter_type())
    , operator_(opertr)
    , inverse_hessian_operator_(inverse_hessian_operator)
  {}

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return operator_.source_space();
  }

  bool linear() const override final
  {
    return false;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_range(range_vector);
    VectorType u_update = range_vector.copy();
    u_update.set_all(0.);
    operator_.apply(source_function, u_update, param);
    inverse_hessian_operator_.apply_inverse_hessian(u_update, range_vector, param);
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return operator_.range_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return operator_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    VectorType u_update = range_vector.copy();
    u_update.set_all(0.);
    operator_.apply(source_vector, u_update, param);
    inverse_hessian_operator_.apply_inverse_hessian(u_update, range_vector, param);
  }

  /// \}

  const OperatorType& operator_;
  const InverseHessianOperatorType& inverse_hessian_operator_;
}; // class EntropicCoordinatesOperator<...>


template <class DensityOperatorImp, class AdvectionOperatorImp, class RhsOperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesCombinedOperator
  : public OperatorInterface<typename AdvectionOperatorImp::AGV,
                             AdvectionOperatorImp::s_r,
                             AdvectionOperatorImp::s_rC,
                             AdvectionOperatorImp::r_r,
                             AdvectionOperatorImp::r_rC,
                             typename AdvectionOperatorImp::F,
                             typename AdvectionOperatorImp::M,
                             typename AdvectionOperatorImp::SGV,
                             typename AdvectionOperatorImp::RGV>
{
  using ThisType = EntropicCoordinatesCombinedOperator;
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::AGV,
                                     AdvectionOperatorImp::s_r,
                                     AdvectionOperatorImp::s_rC,
                                     AdvectionOperatorImp::r_r,
                                     AdvectionOperatorImp::r_rC,
                                     typename AdvectionOperatorImp::F,
                                     typename AdvectionOperatorImp::M,
                                     typename AdvectionOperatorImp::SGV,
                                     typename AdvectionOperatorImp::RGV>;

  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using DensityOperatorType = DensityOperatorImp;
  using AdvectionOperatorType = AdvectionOperatorImp;
  using RhsOperatorType = RhsOperatorImp;
  using InverseHessianOperatorType = InverseHessianOperatorImp;

  EntropicCoordinatesCombinedOperator(const DensityOperatorType& density_op,
                                      const AdvectionOperatorType& advection_op,
                                      const RhsOperatorType& rhs_op,
                                      const InverseHessianOperatorType& inverse_hessian_operator)
    : BaseType(density_op.parameter_type() + advection_op.parameter_type() + rhs_op.parameter_type()
               + inverse_hessian_operator.parameter_type())
    , density_op_(density_op)
    , advection_op_(advection_op)
    , rhs_op_(rhs_op)
    , inverse_hessian_operator_(inverse_hessian_operator)
    , reg_indicators_(advection_op_.source_space().grid_view().size(0), false)
    , u_update_(advection_op_.source_space().mapper().size())
    , rhs_update_(u_update_.size())
  {}

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return advection_op_.range_space();
  }

  bool linear() const override final
  {
    return false;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_range(range_vector);
    density_op_.apply(source_function, range_vector, param);
    u_update_ = range_vector;
    rhs_update_ = range_vector;
    advection_op_.apply(range_vector, u_update_, param);
    u_update_ *= -1.;
    rhs_op_.apply(range_vector, rhs_update_, param);
    u_update_ += rhs_update_;
    std::fill(reg_indicators_.begin(), reg_indicators_.end(), false);
    inverse_hessian_operator_.apply_inverse_hessian(u_update_, reg_indicators_, range_vector, param);
  } // ... apply(...)

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return advection_op_.source_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return advection_op_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    density_op_.apply(source_vector, range_vector, param);
    u_update_ = range_vector;
    rhs_update_ = range_vector;
    advection_op_.apply(range_vector, u_update_, param);
    u_update_ *= -1.;
    rhs_op_.apply(range_vector, rhs_update_, param);
    u_update_ += rhs_update_;
    std::fill(reg_indicators_.begin(), reg_indicators_.end(), false);
    inverse_hessian_operator_.apply_inverse_hessian(u_update_, reg_indicators_, range_vector, param);
  } // ... apply(...)

  /// \}

  const std::vector<bool>& reg_indicators() const
  {
    return reg_indicators_;
  }

  const DensityOperatorType& density_op_;
  const AdvectionOperatorType& advection_op_;
  const RhsOperatorType& rhs_op_;
  const InverseHessianOperatorType& inverse_hessian_operator_;
  mutable std::vector<bool> reg_indicators_;
  mutable VectorType u_update_;
  mutable VectorType rhs_update_;
}; // class EntropicCoordinatesCombinedOperator


template <class AdvectionOperatorType, class ProblemType>
class EntropicCoordinatesMasslumpedOperator
  : public OperatorInterface<typename AdvectionOperatorType::AGV,
                             AdvectionOperatorType::s_r,
                             AdvectionOperatorType::s_rC,
                             AdvectionOperatorType::r_r,
                             AdvectionOperatorType::r_rC,
                             typename AdvectionOperatorType::F,
                             typename AdvectionOperatorType::M,
                             typename AdvectionOperatorType::SGV,
                             typename AdvectionOperatorType::RGV>
{
public:
  using ThisType = EntropicCoordinatesMasslumpedOperator;
  using BaseType = OperatorInterface<typename AdvectionOperatorType::AGV,
                                     AdvectionOperatorType::s_r,
                                     AdvectionOperatorType::s_rC,
                                     AdvectionOperatorType::r_r,
                                     AdvectionOperatorType::r_rC,
                                     typename AdvectionOperatorType::F,
                                     typename AdvectionOperatorType::M,
                                     typename AdvectionOperatorType::SGV,
                                     typename AdvectionOperatorType::RGV>;


  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;
  static constexpr size_t dimDomain = AdvectionOperatorType::SGV::dimension;
  static constexpr size_t dimRange = AdvectionOperatorType::s_r;
  static constexpr double interval_length = 2. / (dimRange - 1.);
  static constexpr size_t first_positive_index_1d = dimRange / 2;
  using RangeFieldType = typename BaseType::FieldType;
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;
  using DomainType = FieldVector<RangeFieldType, dimDomain>;
  using BoundaryFluxesMapType = std::map<DomainType, DynamicVector<RangeFieldType>, XT::Common::FieldVectorFloatLess>;
  using IntersectionType = XT::Grid::extract_intersection_t<typename AdvectionOperatorType::SGV>;

  EntropicCoordinatesMasslumpedOperator(const AdvectionOperatorType& advection_op,
                                        const ProblemType& problem,
                                        const RangeFieldType dx,
                                        const BoundaryFluxesMapType& boundary_fluxes)
    : BaseType(advection_op.parameter_type())
    , advection_op_(advection_op)
    , u_iso_(problem.basis_functions().u_iso())
    , basis_integrated_(problem.basis_functions().integrated())
    , sigma_a_(problem.sigma_a())
    , sigma_s_(problem.sigma_s())
    , Q_(problem.Q())
    , exp_evaluations_(source_space().mapper().size())
    , dummy_range_(dimRange)
    , dx_(dx)
    , h_inv_(1. / dx_)
    , quad_points_(dimDomain * dimRange)
    , quad_weights_(dimRange, 0.)
    , boundary_fluxes_(boundary_fluxes)
    , lower_left_(XT::Common::from_string<DomainType>(problem.grid_config()["lower_left"]))
    , upper_right_(XT::Common::from_string<DomainType>(problem.grid_config()["upper_right"]))
  {
    if constexpr (dimDomain == 1) {
      const auto& nodes = problem.basis_functions().partitioning();
      for (size_t jj = 0; jj < dimRange; ++jj)
        quad_points_[jj] = nodes[jj];
    } else if constexpr (dimDomain == 3) {
      const auto& triangulation = problem.basis_functions().triangulation();
      const auto vertex_quadrature = problem.basis_functions().vertex_quadrature();
      const auto quadratures = triangulation.quadrature_rules(0, vertex_quadrature);
      const auto& vertices = triangulation.vertices();
      std::vector<FieldVector<RangeFieldType, dimDomain>> quad_points(dimRange);
      for (const auto& vertex : vertices) {
        quad_points[vertex->index()] = vertex->position();
        for (size_t dd = 0; dd < dimDomain; ++dd)
          quad_points_[dd * dimRange + vertex->index()] = vertex->position()[dd];
      }
      for (const auto& quadrature : quadratures) {
        for (const auto& point : quadrature) {
          for (size_t jj = 0; jj < dimRange; ++jj) {
            if (XT::Common::FloatCmp::eq(quad_points[jj], point.position())) {
              quad_weights_[jj] += point.weight();
              break;
            }
            if (jj == dimRange - 1)
              DUNE_THROW(Dune::MathError, "Point is not on vertex!");
          } // jj
        } // points
      } // quadratures
    } else {
      DUNE_THROW(Dune::NotImplemented, "Only implemented for dimDomain == 1 or dimDomain == 3");
    }
  } // EntropicCoordinatesMasslumpedOperator(...)

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return advection_op_.range_space();
  }

  bool linear() const override final
  {
    return false;
  }

  // Drop this override to allow for the default implementation in OperatorInterface, which interpolates source_function
  // and calls the other apply below.
  void apply(SourceFunctionType /*source_function*/,
             VectorType& /*range_vector*/,
             const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW(Exceptions::operator_error,
               "Not available for non-discrete functions, call apply(source_vector, range_vector, param) instead!");
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return advection_op_.source_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return advection_operator_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    // Store evaluations of exp(alpha_i)
    assert(source_vector.size() < size_t(std::numeric_limits<int>::max()));
    XT::Common::Mkl::exp(static_cast<int>(source_vector.size()), source_vector.data(), exp_evaluations_.data());
    auto* range_data = range_vector.data();
    const auto& grid_view = source_space().grid_view();
    const size_t num_entities = grid_view.size(0);
    if constexpr (dimDomain == 1) {
      // left boundary flux
      const auto& left_boundary_flux = boundary_fluxes_.at(lower_left_);
      for (size_t jj = 0; jj < dimRange; ++jj)
        range_data[jj] += left_boundary_flux[jj] * h_inv_;
      // inner fluxes and rhs
      for (size_t ii = 0; ii < num_entities; ++ii) {
        const auto offset = ii * dimRange;
        auto* range_entity = range_data + offset;
        auto* range_left = (ii == 0 ? dummy_range_.data() : range_entity - dimRange);
        auto* range_right = (ii == num_entities - 1 ? dummy_range_.data() : range_entity + dimRange);
        const auto* psi = exp_evaluations_.data() + offset;
        // flux
        flux_1d(range_entity, range_left, range_right, psi);
        // rhs
        entity_center_[0] = lower_left_[0] + (ii + 0.5) * dx_;
        rhs_1d(range_entity, psi);
      } // entities
      // left boundary flux
      const auto& right_boundary_flux = boundary_fluxes_.at(upper_right_);
      for (size_t jj = 0; jj < dimRange; ++jj)
        range_data[(num_entities - 1) * dimRange + jj] += right_boundary_flux[jj] * h_inv_;
    } else if constexpr (dimDomain == 3) {
      for (auto&& entity : Dune::elements(grid_view)) {
        const auto entity_index = grid_view.indexSet().index(entity);
        const auto offset = entity_index * dimRange;
        const auto* psi = exp_evaluations_.data() + offset;
        auto* range_entity = range_data + offset;
        // flux
        for (auto&& intersection : Dune::intersections(grid_view, entity)) {
          const bool boundary = !intersection.neighbor();
          const auto& outside_entity = intersection.outside();
          auto* range_outside =
              boundary ? dummy_range_.data() : range_data + grid_view.indexSet().index(outside_entity) * dimRange;
          flux_3d(range_entity, range_outside, psi, intersection);
        } // intersections
        // rhs
        entity_center_ = entity.geometry().center();
        rhs_3d(range_entity, psi);
      } // entities
    }
    // inverse Hessian
    const auto* psi = exp_evaluations_.data();
    if constexpr (dimDomain == 1)
      apply_inverse_hessian_1d(num_entities, range_data, psi);
    else
      apply_inverse_hessian_3d(num_entities, range_data, psi);
  } // void apply(...)

  /// \}

  void flux_1d(RangeFieldType* BOOST_RESTRICT range_entity,
               RangeFieldType* BOOST_RESTRICT range_left,
               RangeFieldType* BOOST_RESTRICT range_right,
               const RangeFieldType* BOOST_RESTRICT psi) const
  {
    // fluxes in negative direction
    for (size_t jj = 0; jj < first_positive_index_1d; ++jj) {
      const auto weight = (jj == 0 ? interval_length * 0.5 : interval_length);
      const auto flux_jj = -quad_points_[jj] * psi[jj] * h_inv_ * weight;
      range_entity[jj] -= flux_jj;
      range_left[jj] += flux_jj;
    }
    // fluxes in positive direction
    for (size_t jj = first_positive_index_1d; jj < dimRange; ++jj) {
      const auto weight = (jj == dimRange - 1 ? interval_length * 0.5 : interval_length);
      const auto flux_jj = quad_points_[jj] * psi[jj] * h_inv_ * weight;
      range_entity[jj] -= flux_jj;
      range_right[jj] += flux_jj;
    }
  }

  void flux_3d(RangeFieldType* BOOST_RESTRICT range_entity,
               RangeFieldType* BOOST_RESTRICT range_outside,
               const RangeFieldType* BOOST_RESTRICT psi,
               const IntersectionType& intersection) const
  {
    const bool boundary = !intersection.neighbor();
    const auto dd = intersection.indexInInside() / 2;
    const auto n_dd = intersection.centerUnitOuterNormal()[dd];
    const auto quad_offset = dd * dimRange;
    for (size_t jj = 0; jj < dimRange; ++jj) {
      const auto v_times_n = quad_points_[quad_offset + jj] * n_dd;
      if (v_times_n > 0.) {
        const auto flux_jj = psi[jj] * quad_weights_[jj] * h_inv_ * v_times_n;
        range_entity[jj] -= flux_jj;
        range_outside[jj] += flux_jj;
      }
    } // jj
    if (boundary) {
      const auto& boundary_flux = boundary_fluxes_.at(intersection.geometry().center());
      for (size_t jj = 0; jj < dimRange; ++jj)
        range_entity[jj] += boundary_flux[jj] * h_inv_;
    }
  }

  // Note: entity_center_ has to be set to entity.geometry().center() before calling this function
  void rhs_1d(RangeFieldType* BOOST_RESTRICT range_entity, const RangeFieldType* BOOST_RESTRICT psi) const
  {
    const auto sigma_a_value = sigma_a_->evaluate(entity_center_)[0];
    const auto sigma_s_value = sigma_s_->evaluate(entity_center_)[0];
    const auto Q_value = Q_->evaluate(entity_center_)[0];
    auto density = XT::Common::reduce(psi + 1, psi + dimRange - 1, 0.);
    density += 0.5 * (psi[0] + psi[dimRange - 1]);
    density *= interval_length;
    const auto sigma_s_factor = sigma_s_value * density;
    const auto sigma_t_factor = (sigma_a_value + sigma_s_value) * interval_length;
    range_entity[0] += u_iso_[0] * sigma_s_factor + basis_integrated_[0] * Q_value - psi[0] * sigma_t_factor * 0.5;
    for (size_t jj = 1; jj < dimRange - 1; ++jj)
      range_entity[jj] += u_iso_[jj] * sigma_s_factor + basis_integrated_[jj] * Q_value - psi[jj] * sigma_t_factor;
    range_entity[dimRange - 1] += u_iso_[dimRange - 1] * sigma_s_factor + basis_integrated_[dimRange - 1] * Q_value
                                  - psi[dimRange - 1] * sigma_t_factor * 0.5;
  }

  // Note: entity_center_ has to be set to entity.geometry().center() before calling this function
  void rhs_3d(RangeFieldType* BOOST_RESTRICT range_entity, const RangeFieldType* BOOST_RESTRICT psi) const
  {
    const auto sigma_a_value = sigma_a_->evaluate(entity_center_)[0];
    const auto sigma_s_value = sigma_s_->evaluate(entity_center_)[0];
    const auto sigma_t_value = sigma_a_value + sigma_s_value;
    const auto Q_value = Q_->evaluate(entity_center_)[0];
    const auto density = XT::Common::transform_reduce(psi, psi + dimRange, quad_weights_.begin(), 0.);
    const auto sigma_s_factor = sigma_s_value * density;
    for (size_t jj = 0; jj < dimRange; ++jj)
      range_entity[jj] +=
          u_iso_[jj] * sigma_s_factor + basis_integrated_[jj] * Q_value - psi[jj] * quad_weights_[jj] * sigma_t_value;
  }

  static void apply_inverse_hessian_1d(const size_t num_entities,
                                       RangeFieldType* BOOST_RESTRICT range_data,
                                       const RangeFieldType* BOOST_RESTRICT psi)
  {
    for (size_t ii = 0; ii < num_entities; ++ii) {
      const auto offset = ii * dimRange;
      range_data[offset] /= psi[offset] * interval_length * 0.5;
      for (size_t jj = 1; jj < dimRange - 1; ++jj)
        range_data[offset + jj] /= psi[offset + jj] * interval_length;
      range_data[offset + dimRange - 1] /= psi[offset + dimRange - 1] * interval_length * 0.5;
    } // ii
  }

  void apply_inverse_hessian_3d(const size_t num_entities,
                                RangeFieldType* BOOST_RESTRICT range_data,
                                const RangeFieldType* BOOST_RESTRICT psi) const
  {
    for (size_t ii = 0; ii < num_entities; ++ii) {
      const auto offset = ii * dimRange;
      for (size_t jj = 0; jj < dimRange; ++jj)
        range_data[offset + jj] /= psi[offset + jj] * quad_weights_[jj];
    } // ii
  }

  const AdvectionOperatorType& advection_op_;
  DynamicVector<RangeFieldType> u_iso_;
  DynamicVector<RangeFieldType> basis_integrated_;
  std::unique_ptr<ParameterFunctionType> sigma_a_;
  std::unique_ptr<ParameterFunctionType> sigma_s_;
  std::unique_ptr<ParameterFunctionType> Q_;
  mutable std::vector<RangeFieldType> exp_evaluations_;
  mutable DomainType entity_center_;
  mutable std::vector<RangeFieldType> dummy_range_;
  const RangeFieldType dx_;
  const RangeFieldType h_inv_;
  std::vector<RangeFieldType> quad_points_;
  std::vector<RangeFieldType> quad_weights_;
  const BoundaryFluxesMapType& boundary_fluxes_;
  const DomainType lower_left_;
  const DomainType upper_right_;
}; // class EntropicCoordinatesMasslumpedOperator


template <class AdvectionOperatorImp, class EntropySolverImp>
class EntropyBasedMomentFvOperator
  : public OperatorInterface<typename AdvectionOperatorImp::AGV,
                             AdvectionOperatorImp::s_r,
                             AdvectionOperatorImp::s_rC,
                             AdvectionOperatorImp::r_r,
                             AdvectionOperatorImp::r_rC,
                             typename AdvectionOperatorImp::F,
                             typename AdvectionOperatorImp::M,
                             typename AdvectionOperatorImp::SGV,
                             typename AdvectionOperatorImp::RGV>
{
public:
  using ThisType = EntropyBasedMomentFvOperator;
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::AGV,
                                     AdvectionOperatorImp::s_r,
                                     AdvectionOperatorImp::s_rC,
                                     AdvectionOperatorImp::r_r,
                                     AdvectionOperatorImp::r_rC,
                                     typename AdvectionOperatorImp::F,
                                     typename AdvectionOperatorImp::M,
                                     typename AdvectionOperatorImp::SGV,
                                     typename AdvectionOperatorImp::RGV>;

  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AdvectionOperatorType = AdvectionOperatorImp;
  using EntropySolverType = EntropySolverImp;

  EntropyBasedMomentFvOperator(const AdvectionOperatorType& advection_operator, const EntropySolverType& entropy_solver)
    : BaseType(advection_operator.parameter_type())
    , advection_operator_(advection_operator)
    , entropy_solver_(entropy_solver)
  {}

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  bool linear() const override final
  {
    return false;
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  // Drop this override to allow for the default implementation in OperatorInterface, which interpolates source_function
  // and calls the other apply below. If entropy_solver_ supports apply(source_function, range_vector), simply fill this
  // apply() with the code from the other one
  void apply(SourceFunctionType /*source_function*/,
             VectorType& /*range_vector*/,
             const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW(Exceptions::operator_error,
               "Not available for non-discrete functions, call apply(source_vector, range_vector, param) instead!");
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return advection_operator_.source_space();
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return advection_operator_.assembly_grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    // solve optimization problems and regularize if necessary
    VectorType regularized = range_vector;
    entropy_solver_.apply(source_vector, regularized, param);

    range_vector.set_all(0.);
    advection_operator_.apply(regularized, range_vector, param);
  } // ... apply(...)

  /// \}

  const AdvectionOperatorType& advection_operator_;
  const EntropySolverType& entropy_solver_;
}; // class EntropyBasedMomentFvOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
