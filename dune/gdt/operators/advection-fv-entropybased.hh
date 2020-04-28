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

#include <dune/xt/common/numeric.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {


template <class OperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesOperator
  : public OperatorInterface<typename OperatorImp::MatrixType, typename OperatorImp::SGV, OperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename OperatorImp::MatrixType, typename OperatorImp::SGV, OperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using OperatorType = OperatorImp;
  using InverseHessianOperatorType = InverseHessianOperatorImp;

  EntropicCoordinatesOperator(const OperatorType& operator_in,
                              const InverseHessianOperatorType& inverse_hessian_operator)
    : operator_(operator_in)
    , inverse_hessian_operator_(inverse_hessian_operator)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return operator_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return operator_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    VectorType u_update = range;
    std::fill(u_update.begin(), u_update.end(), 0.);
    operator_.apply(source, u_update, param);
    inverse_hessian_operator_.apply_inverse_hessian(u_update, range, param);
  }

  const OperatorType& operator_;
  const InverseHessianOperatorType& inverse_hessian_operator_;
}; // class EntropicCoordinatesOperator<...>

template <class DensityOperatorImp, class AdvectionOperatorImp, class RhsOperatorImp, class InverseHessianOperatorImp>
class EntropicCoordinatesCombinedOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SGV,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SGV,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
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
    : density_op_(density_op)
    , advection_op_(advection_op)
    , rhs_op_(rhs_op)
    , inverse_hessian_operator_(inverse_hessian_operator)
    , reg_indicators_(advection_op_.source_space().grid_view().size(0), false)
    , u_update_(advection_op_.source_space().mapper().size())
    , rhs_update_(u_update_.size())
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return advection_op_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_op_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    density_op_.apply(source, range, param);
    u_update_ = range;
    rhs_update_ = range;
    advection_op_.apply(range, u_update_, param);
    u_update_ *= -1.;
    rhs_op_.apply(range, rhs_update_, param);
    u_update_ += rhs_update_;
    std::fill(reg_indicators_.begin(), reg_indicators_.end(), false);
    inverse_hessian_operator_.apply_inverse_hessian(u_update_, reg_indicators_, range, param);
  }

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
}; // class EntropicCoordinatesCombinedOperator<...>


template <class AdvectionOperatorType, class ProblemType>
class EntropicCoordinatesMasslumpedOperator
  : public OperatorInterface<typename AdvectionOperatorType::MatrixType,
                             typename AdvectionOperatorType::SGV,
                             AdvectionOperatorType::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorType::MatrixType,
                                     typename AdvectionOperatorType::SGV,
                                     AdvectionOperatorType::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;
  static constexpr size_t dimDomain = AdvectionOperatorType::SGV::dimension;
  static constexpr size_t dimRange = AdvectionOperatorType::s_r;
  static constexpr double interval_length = 2. / (dimRange - 1.);
  using RangeFieldType = typename BaseType::FieldType;
  using ParameterFunctionType = typename ProblemType::ScalarFunctionType;
  using DomainType = FieldVector<RangeFieldType, dimDomain>;
  using BoundaryFluxesMapType = std::map<DomainType, DynamicVector<RangeFieldType>, XT::Common::FieldVectorFloatLess>;

  EntropicCoordinatesMasslumpedOperator(const AdvectionOperatorType& advection_op,
                                        const ProblemType& problem,
                                        const RangeFieldType dx,
                                        const BoundaryFluxesMapType& boundary_fluxes)
    : advection_op_(advection_op)
    , u_iso_(problem.basis_functions().u_iso())
    , basis_integrated_(problem.basis_functions().integrated())
    , sigma_a_(problem.sigma_a())
    , sigma_s_(problem.sigma_s())
    , Q_(problem.Q())
    , exp_evaluations_(source_space().mapper().size())
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
  }

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return advection_op_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_op_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/) const override final
  {
    // Store evaluations of exp(alpha_i)
    assert(source.size() < std::numeric_limits<int>::max());
    XT::Common::Mkl::exp(static_cast<int>(source.size()), source.data(), exp_evaluations_.data());
    auto* range_data = range.data();
    std::fill(range_data, range_data + range.size(), 0.);
    const auto& grid_view = source_space().grid_view();
    const size_t num_entities = grid_view.size(0);
    DomainType entity_center;
    if constexpr (dimDomain == 1) {
      for (size_t ii = 0; ii < num_entities; ++ii) {
        const bool left_boundary = (ii == 0);
        const bool right_boundary = (ii == num_entities - 1);
        const auto offset = ii * dimRange;
        auto* range_entity = range_data + offset;
        auto* range_left = range_entity - dimRange;
        auto* range_right = range_entity + dimRange;
        const auto* psi = exp_evaluations_.data() + offset;
        const auto* boundary_flux =
            left_boundary ? &(boundary_fluxes_.at(lower_left_)) : &(boundary_fluxes_.at(upper_right_));
        entity_center[0] = lower_left_[0] + (ii + 0.5) * dx_;
        const auto sigma_a_value = sigma_a_->evaluate(entity_center)[0];
        const auto sigma_s_value = sigma_s_->evaluate(entity_center)[0];
        const auto sigma_t_value = sigma_a_value + sigma_s_value;
        const auto Q_value = Q_->evaluate(entity_center)[0];
        auto density = XT::Common::reduce(psi + 1, psi + dimRange - 1, 0.);
        density += 0.5 * (psi[0] + psi[dimRange - 1]);
        density *= interval_length;
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // flux
          const auto& v = quad_points_[jj];
          const auto weight = (jj == 0 || jj == dimRange - 1) ? interval_length / 2. : interval_length;
          const auto flux_jj = std::abs(v) * psi[jj] * h_inv_ * weight;
          range_entity[jj] -= flux_jj;
          if (v > 0.) {
            if (!right_boundary)
              range_right[jj] += flux_jj;
          } else if (!left_boundary) {
            range_left[jj] += flux_jj;
          }
          if (left_boundary || right_boundary)
            range_entity[jj] += (*boundary_flux)[jj] * h_inv_;
          // rhs
          range_entity[jj] +=
              sigma_s_value * density * u_iso_[jj] + Q_value * basis_integrated_[jj] - psi[jj] * weight * sigma_t_value;
        } // jj
      } // entities
    } else if constexpr (dimDomain == 3) {
      for (auto&& entity : Dune::elements(grid_view)) {
        const auto entity_index = grid_view.indexSet().index(entity);
        const auto offset = entity_index * dimRange;
        const auto* psi = exp_evaluations_.data() + offset;
        auto* range_entity = range_data + offset;
        // flux
        for (auto&& intersection : Dune::intersections(grid_view, entity)) {
          const bool boundary = !intersection.neighbor();
          const auto& outside_entity = boundary ? entity : intersection.outside();
          const auto outside_index = grid_view.indexSet().index(outside_entity);
          auto* range_outside = range_data + outside_index * dimRange;
          const auto dd = intersection.indexInInside() / 2;
          const auto& n = intersection.centerUnitOuterNormal();
          const auto* boundary_flux = boundary ? &(boundary_fluxes_.at(intersection.geometry().center())) : nullptr;
          for (size_t jj = 0; jj < dimRange; ++jj) {
            const auto v_times_n = quad_points_[dd * dimRange + jj] * n[dd];
            if (v_times_n > 0.) {
              const auto flux_jj = v_times_n * psi[jj] * quad_weights_[jj] * h_inv_;
              range_entity[jj] -= flux_jj;
              if (!boundary)
                range_outside[jj] += flux_jj;
            } else if (boundary) {
              range_entity[jj] += (*boundary_flux)[jj] * h_inv_;
            }
          } // jj
        } // intersections
        // rhs
        entity_center = entity.geometry().center();
        const auto sigma_a_value = sigma_a_->evaluate(entity_center)[0];
        const auto sigma_s_value = sigma_s_->evaluate(entity_center)[0];
        const auto sigma_t_value = sigma_a_value + sigma_s_value;
        const auto Q_value = Q_->evaluate(entity_center)[0];
        const auto density = XT::Common::transform_reduce(psi, psi + dimRange, quad_weights_.begin(), 0.);
        for (size_t jj = 0; jj < dimRange; ++jj) {
          range_entity[jj] += sigma_s_value * density * u_iso_[jj] + Q_value * basis_integrated_[jj]
                              - psi[jj] * quad_weights_[jj] * sigma_t_value;
        } // jj
      } // entities
    }
    // inverse Hessian
    const auto* psi = exp_evaluations_.data();
    for (size_t ii = 0; ii < num_entities; ++ii) {
      const auto offset = ii * dimRange;
      for (size_t jj = 0; jj < dimRange; ++jj) {
        auto& val = range_data[offset + jj];
        if constexpr (dimDomain == 1) {
          const auto weight = (jj == 0 || jj == dimRange - 1) ? interval_length / 2. : interval_length;
          val /= psi[offset + jj] * weight;
        } else {
          val /= psi[offset + jj] * quad_weights_[jj];
        }
        if (std::isnan(val) || std::isinf(val))
          DUNE_THROW(Dune::MathError, "inf or nan in range!");
      } // jj
    } // ii
  } // void apply(...)

  const AdvectionOperatorType& advection_op_;
  DynamicVector<RangeFieldType> u_iso_;
  DynamicVector<RangeFieldType> basis_integrated_;
  std::unique_ptr<ParameterFunctionType> sigma_a_;
  std::unique_ptr<ParameterFunctionType> sigma_s_;
  std::unique_ptr<ParameterFunctionType> Q_;
  mutable std::vector<RangeFieldType> exp_evaluations_;
  const RangeFieldType dx_;
  const RangeFieldType h_inv_;
  std::vector<RangeFieldType> quad_points_;
  std::vector<RangeFieldType> quad_weights_;
  const BoundaryFluxesMapType& boundary_fluxes_;
  const DomainType lower_left_;
  const DomainType upper_right_;
}; // class EntropicCoordinatesMasslumpedOperator<...>


template <class AdvectionOperatorImp, class EntropySolverImp>
class EntropyBasedMomentFvOperator
  : public OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                             typename AdvectionOperatorImp::SGV,
                             AdvectionOperatorImp::s_r>
{
  using BaseType = OperatorInterface<typename AdvectionOperatorImp::MatrixType,
                                     typename AdvectionOperatorImp::SGV,
                                     AdvectionOperatorImp::s_r>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using AdvectionOperatorType = AdvectionOperatorImp;
  using EntropySolverType = EntropySolverImp;

  EntropyBasedMomentFvOperator(const AdvectionOperatorType& advection_operator, const EntropySolverType& entropy_solver)
    : advection_operator_(advection_operator)
    , entropy_solver_(entropy_solver)
  {}

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return advection_operator_.source_space();
  }

  const RangeSpaceType& range_space() const override final
  {
    return advection_operator_.range_space();
  }

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param) const override final
  {
    // solve optimization problems and regularize if necessary
    VectorType regularized = range;
    entropy_solver_.apply(source, regularized, param);

    std::fill(range.begin(), range.end(), 0.);
    advection_operator_.apply(regularized, range, param);
  }

  const AdvectionOperatorType& advection_operator_;
  const EntropySolverType& entropy_solver_;
}; // class EntropyBasedMomentFvOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_ENTROPYBASED_HH
