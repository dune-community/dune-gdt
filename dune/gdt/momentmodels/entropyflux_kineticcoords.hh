// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_KINETICCOORDS_HH
#define DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_KINETICCOORDS_HH

#include <list>
#include <memory>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/functions/interfaces/flux-function.hh>

#include <dune/gdt/momentmodels/basisfunctions.hh>
#include <dune/gdt/momentmodels/entropyflux_implementations.hh>
#include <dune/gdt/momentmodels/entropyflux.hh>

namespace Dune {
namespace GDT {


template <class GridViewImp, class MomentBasisImp>
class EntropyBasedFluxEntropyCoordsFunction
  : public XT::Functions::FluxFunctionInterface<XT::Grid::extract_entity_t<GridViewImp>,
                                                MomentBasisImp::dimRange,
                                                MomentBasisImp::dimFlux,
                                                MomentBasisImp::dimRange,
                                                typename MomentBasisImp::R>
{
  using BaseType = typename XT::Functions::FluxFunctionInterface<XT::Grid::extract_entity_t<GridViewImp>,
                                                                 MomentBasisImp::dimRange,
                                                                 MomentBasisImp::dimFlux,
                                                                 MomentBasisImp::dimRange,
                                                                 typename MomentBasisImp::R>;
  using ThisType = EntropyBasedFluxEntropyCoordsFunction;

public:
  using GridViewType = GridViewImp;
  using MomentBasis = MomentBasisImp;
  using IndexSetType = typename GridViewType::IndexSet;
  static const size_t dimFlux = MomentBasis::dimFlux;
  static const size_t basis_dimRange = MomentBasis::dimRange;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::LocalFunctionType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::StateType;
  using ImplementationType = EntropyBasedFluxImplementation<MomentBasis>;
  using AlphaReturnType = typename ImplementationType::AlphaReturnType;
  using VectorType = typename ImplementationType::VectorType;
  using I = XT::Grid::extract_intersection_t<GridViewType>;

  explicit EntropyBasedFluxEntropyCoordsFunction(
      const MomentBasis& basis_functions,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52))
    : implementation_(std::make_shared<ImplementationType>(
          basis_functions, tau, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon))
  {}

  explicit EntropyBasedFluxEntropyCoordsFunction(EntropyBasedFluxFunction<GridViewType, MomentBasis>& other)
    : implementation_(other.implementation_)
  {}


  static const constexpr bool available = true;

  class Localfunction : public LocalFunctionType
  {
    using BaseType = LocalFunctionType;

  public:
    using typename BaseType::DynamicJacobianRangeType;
    using typename BaseType::E;
    using typename BaseType::RangeReturnType;

    Localfunction(const ImplementationType& implementation)
      : implementation_(implementation)
    {}

    virtual int order(const XT::Common::Parameter&) const override final
    {
      return 1.;
    }

    virtual RangeReturnType evaluate(const DomainType& /*point_in_reference_element*/,
                                     const StateType& alpha,
                                     const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      return implementation_.evaluate_with_alpha(alpha);
    }

    virtual void jacobian(const DomainType& /*point_in_reference_element*/,
                          const StateType& alpha,
                          DynamicJacobianRangeType& result,
                          const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      implementation_.jacobian_with_alpha(alpha, result);
    } // ... jacobian(...)

  private:
    const ImplementationType& implementation_;
  }; // class Localfunction

  virtual bool x_dependent() const override final
  {
    return false;
  }

  virtual std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<Localfunction>(*implementation_);
  }

  virtual std::unique_ptr<Localfunction> derived_local_function() const
  {
    return std::make_unique<Localfunction>(*implementation_);
  }

  StateType evaluate_kinetic_flux_precomputed(const StateType& flux_1,
                                              const StateType& flux_2,
                                              const DomainType& n_ij,
                                              const size_t dd) const
  {
    return (flux_1 + flux_2) * n_ij[dd];
  } // StateType evaluate_kinetic_flux(...)

  StateType evaluate_kinetic_flux(const E& /*inside_entity*/,
                                  const E& /*outside_entity*/,
                                  const StateType& alpha_i,
                                  const StateType& alpha_j,
                                  const DomainType& n_ij,
                                  const size_t dd) const
  {
    return implementation_->evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // StateType evaluate_kinetic_flux(...)

  const MomentBasis& basis_functions() const
  {
    return implementation_->basis_functions();
  }

  StateType get_u(const StateType& alpha) const
  {
    return implementation_->get_u(alpha);
  }

  StateType get_alpha(const StateType& u) const
  {
    const auto alpha = implementation_->get_alpha(u)->first;
    StateType ret;
    std::copy(alpha.begin(), alpha.end(), ret.begin());
    return ret;
  }

  template <class FluxesMapType>
  void calculate_minmod_density_reconstruction(const FieldVector<StateType, 3>& alphas,
                                               FluxesMapType& precomputed_fluxes,
                                               const size_t dd) const
  {
    implementation_->calculate_minmod_density_reconstruction(alphas, precomputed_fluxes, dd);
  }

  StateType calculate_boundary_flux(const StateType& alpha, const I& intersection)
  {
    return implementation_->calculate_boundary_flux(alpha, intersection);
  }

  void apply_inverse_hessian(const StateType& alpha, const StateType& u, StateType& Hinv_u) const
  {
    implementation_->apply_inverse_hessian(alpha, u, Hinv_u);
  }

private:
  std::shared_ptr<ImplementationType> implementation_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_KINETICCOORDS_HH
