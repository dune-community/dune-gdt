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

#include <dune/xt/functions/checkerboard.hh>

#include <dune/gdt/assembler/system.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType, class BasisfunctionType>
class KineticIsotropicLocalFunctor
  : public XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  using BaseType = typename XT::Grid::Functor::Codim0<typename DiscreteFunctionType::SpaceType::GridLayerType>;
  using CheckerboardType = XT::Functions::CheckerboardFunction<typename DiscreteFunctionType::EntityType,
                                                               typename DiscreteFunctionType::DomainFieldType,
                                                               DiscreteFunctionType::dimDomain,
                                                               typename DiscreteFunctionType::RangeFieldType,
                                                               1,
                                                               1>;
  using RangeType = typename DiscreteFunctionType::RangeType;
  static constexpr size_t dimRange = DiscreteFunctionType::dimRange;

public:
  using typename BaseType::EntityType;

  KineticIsotropicLocalFunctor(const BasisfunctionType& basis_functions,
                               DiscreteFunctionType& solution,
                               const double dt,
                               const CheckerboardType& sigma_a,
                               const CheckerboardType& sigma_s,
                               const CheckerboardType& Q)

    : basis_functions_(basis_functions)
    , solution_(solution)
    , dt_(dt)
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
    , basis_integrated_(basis_functions_.integrated())
    , u_iso_(basis_functions_.u_iso())
  {}

  virtual void apply_local(const EntityType& entity)
  {
    auto solution_local = solution_.local_discrete_function(entity);

    // get u
    const auto center = entity.geometry().local(entity.geometry().center());
    const auto u0 = solution_local->evaluate(center);
    const auto sigma_a = sigma_a_.local_function(entity)->evaluate(center)[0];
    const auto sigma_s = sigma_s_.local_function(entity)->evaluate(center)[0];
    const auto Q = Q_.local_function(entity)->evaluate(center)[0];
    const auto exp_sigma_a = std::exp(-sigma_a * dt_);
    const auto exp_sigma_s = std::exp(-sigma_s * dt_);
    const auto u_iso = u_iso_ * basis_functions_.density(u0);

    auto u = exp_sigma_a * (u0 * exp_sigma_s + u_iso * (1 - exp_sigma_s));
    if (!XT::Common::is_zero(sigma_a))
      u += basis_integrated_ * (1 - exp_sigma_a) / sigma_a * Q;
    else
      u += basis_integrated_ * dt_ * Q;

    // write to return vector
    auto& local_vector = solution_local->vector();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_vector.set(ii, u[ii]);
  }

private:
  const BasisfunctionType& basis_functions_;
  DiscreteFunctionType& solution_;
  const double dt_;
  const CheckerboardType& sigma_a_;
  const CheckerboardType& sigma_s_;
  const CheckerboardType& Q_;
  const RangeType basis_integrated_;
  const RangeType u_iso_;
};


/** \brief Time stepper solving linear equation d_t u = Au + b by matrix exponential
 */
template <class DiscreteFunctionImp, class BasisfunctionType>
class KineticIsotropicTimeStepper : public TimeStepperInterface<DiscreteFunctionImp>
{
  typedef KineticIsotropicTimeStepper ThisType;
  typedef TimeStepperInterface<DiscreteFunctionImp> BaseType;

public:
  using typename BaseType::DiscreteFunctionType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::EntityType;
  using typename BaseType::RangeFieldType;

  static const size_t dimDomain = DiscreteFunctionType::dimDomain;
  static const size_t dimRange = DiscreteFunctionType::dimRange;
  using CheckerboardType =
      XT::Functions::CheckerboardFunction<EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1>;
  using BaseType::current_solution;
  using BaseType::current_time;

  KineticIsotropicTimeStepper(const BasisfunctionType& basis_functions,
                              DiscreteFunctionType& initial_values,
                              const CheckerboardType& sigma_a,
                              const CheckerboardType& sigma_s,
                              const CheckerboardType& Q,
                              const RangeFieldType t_0 = 0.0)
    : BaseType(t_0, initial_values)
    , basis_functions_(basis_functions)
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
  {}

  virtual RangeFieldType step(const RangeFieldType dt, const RangeFieldType max_dt) override final
  {
    const RangeFieldType actual_dt = std::min(dt, max_dt);
    auto& t = current_time();
    auto& u_n = current_solution();
    KineticIsotropicLocalFunctor<DiscreteFunctionType, BasisfunctionType> functor(
        basis_functions_, u_n, actual_dt, sigma_a_, sigma_s_, Q_);
    SystemAssembler<typename DiscreteFunctionType::SpaceType> assembler(u_n.space());
    assembler.append(functor);
    assembler.assemble(true);
    t += actual_dt;
    return dt;
  } // ... step(...)

private:
  const BasisfunctionType& basis_functions_;
  const CheckerboardType sigma_a_;
  const CheckerboardType sigma_s_;
  const CheckerboardType Q_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TIMESTEPPER_KINETIC_ISOTROPIC_HH
