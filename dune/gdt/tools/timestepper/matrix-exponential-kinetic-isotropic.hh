// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_TIMESTEPPER_KINETIC_ISOTROPIC_HH
#define DUNE_GDT_TIMESTEPPER_KINETIC_ISOTROPIC_HH

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/xt/functions/interfaces/function.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType, class MomentBasis>
class KineticIsotropicLocalFunctor
  : public XT::Grid::ElementFunctor<typename DiscreteFunctionType::SpaceType::GridViewType>
{
  using GridViewType = typename DiscreteFunctionType::SpaceType::GridViewType;
  using BaseType = typename XT::Grid::ElementFunctor<GridViewType>;
  using RangeType = typename DiscreteFunctionType::LocalFunctionType::RangeReturnType;
  using ScalarFunctionType =
      XT::Functions::FunctionInterface<MomentBasis::dimDomain, 1, 1, typename MomentBasis::RangeFieldType>;
  static constexpr size_t dimRange = DiscreteFunctionType::r;

public:
  using typename BaseType::E;

  KineticIsotropicLocalFunctor(const MomentBasis& basis_functions,
                               DiscreteFunctionType& solution,
                               const double dt,
                               const ScalarFunctionType& sigma_a,
                               const ScalarFunctionType& sigma_s,
                               const ScalarFunctionType& Q)

    : basis_functions_(basis_functions)
    , solution_(solution)
    , local_solution_(solution_.local_discrete_function())
    , dt_(dt)
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
    , basis_integrated_(basis_functions_.integrated())
    , u_iso_(basis_functions_.u_iso())
  {}

  KineticIsotropicLocalFunctor(const KineticIsotropicLocalFunctor& other)
    : BaseType(other)
    , basis_functions_(other.basis_functions_)
    , solution_(other.solution_)
    , local_solution_(solution_.local_discrete_function())
    , dt_(other.dt_)
    , sigma_a_(other.sigma_a_)
    , sigma_s_(other.sigma_s_)
    , Q_(other.Q_)
    , basis_integrated_(other.basis_integrated_)
    , u_iso_(other.u_iso_)
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new KineticIsotropicLocalFunctor(*this);
  }

  void apply_local(const E& entity) override final
  {
    local_solution_->bind(entity);

    // get u
    const auto center = entity.geometry().center();
    const auto center_local = entity.geometry().local(center);
    const auto u0 = local_solution_->evaluate(center_local);
    const auto sigma_a = sigma_a_.evaluate(center)[0];
    const auto sigma_s = sigma_s_.evaluate(center)[0];
    const auto Q = Q_.evaluate(center)[0];
    const auto exp_sigma_a = std::exp(-sigma_a * dt_);
    const auto exp_sigma_s = std::exp(-sigma_s * dt_);
    const auto u_iso = u_iso_ * basis_functions_.density(u0);

    auto u = exp_sigma_a * (u0 * exp_sigma_s + u_iso * (1 - exp_sigma_s));
    if (!XT::Common::is_zero(sigma_a))
      u += basis_integrated_ * (1 - exp_sigma_a) / sigma_a * Q;
    else
      u += basis_integrated_ * dt_ * Q;

    // write to return vector
    auto& local_vector = local_solution_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_vector.set_entry(ii, u[ii]);
  }

private:
  const MomentBasis& basis_functions_;
  DiscreteFunctionType& solution_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_solution_;
  const double dt_;
  const ScalarFunctionType& sigma_a_;
  const ScalarFunctionType& sigma_s_;
  const ScalarFunctionType& Q_;
  const RangeType basis_integrated_;
  const RangeType u_iso_;
};


/** \brief Time stepper solving linear equation d_t u = Au + b by matrix exponential
 */
template <class DiscreteFunctionImp, class MomentBasis>
class KineticIsotropicTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  using ThisType = KineticIsotropicTimeStepper;
  using BaseType = TimeStepperInterface<DiscreteFunctionImp>;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeFieldType;
  static constexpr size_t dimDomain = DiscreteFunctionType::d;
  static constexpr size_t dimRange = DiscreteFunctionType::r;
  using ScalarFunctionType = XT::Functions::FunctionInterface<dimDomain, 1, 1, RangeFieldType>;

  using BaseType::current_solution;
  using BaseType::current_time;

  KineticIsotropicTimeStepper(const MomentBasis& basis_functions,
                              DiscreteFunctionType& initial_values,
                              const ScalarFunctionType& sigma_a,
                              const ScalarFunctionType& sigma_s,
                              const ScalarFunctionType& Q,
                              const RangeFieldType t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , basis_functions_(basis_functions)
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
  {}

  RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    KineticIsotropicLocalFunctor<DiscreteFunctionType, MomentBasis> functor(
        basis_functions_, u_n, actual_dt, sigma_a_, sigma_s_, Q_);
    auto walker = XT::Grid::Walker<typename DiscreteFunctionType::SpaceType::GridViewType>(u_n.space().grid_view());
    walker.append(functor);
    walker.walk(true);
    t += actual_dt;
    return dt;
  } // ... step(...)

private:
  const MomentBasis& basis_functions_;
  const ScalarFunctionType& sigma_a_;
  const ScalarFunctionType& sigma_s_;
  const ScalarFunctionType& Q_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_KINETIC_ISOTROPIC_HH
