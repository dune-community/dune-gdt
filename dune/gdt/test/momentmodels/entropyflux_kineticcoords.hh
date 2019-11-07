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

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux_implementations.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>

namespace Dune {
namespace GDT {


template <class GridViewImp, class MomentBasisImp, SlopeLimiterType slope>
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
  using QuadratureWeightsType = typename ImplementationType::QuadratureWeightsType;
  using BoundaryQuadratureWeightsType =
      std::vector<XT::Common::FieldVector<XT::Common::FieldVector<QuadratureWeightsType, 2>, dimFlux>>;
  static const EntropyType entropy = MomentBasis::entropy;

  explicit EntropyBasedFluxEntropyCoordsFunction(
      const MomentBasis& basis_functions,
      const bool disable_realizability_check = false,
      const RangeFieldType tau = 1e-9,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52))
    : implementation_(std::make_shared<ImplementationType>(
          basis_functions, tau, disable_realizability_check, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon))
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

    int order(const XT::Common::Parameter&) const override final
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

  bool x_dependent() const override final
  {
    return false;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<Localfunction>(*implementation_);
  }

  virtual std::unique_ptr<Localfunction> derived_local_function() const
  {
    return std::make_unique<Localfunction>(*implementation_);
  }

  StateType evaluate_kinetic_flux(const E& /*inside_entity*/,
                                  const E& /*outside_entity*/,
                                  const StateType& flux_i,
                                  const StateType& flux_j,
                                  const DomainType& n_ij,
                                  const size_t dd) const
  {
    return (flux_i + flux_j) * n_ij[dd];
  } // StateType evaluate_kinetic_flux(...)

  const MomentBasis& basis_functions() const
  {
    return implementation_->basis_functions();
  }

  StateType get_u(const StateType& alpha) const
  {
    return implementation_->get_u(alpha);
  }

  StateType get_u(const size_t entity_index) const
  {
    return implementation_->get_u((*eta_ast_prime_evaluations_)[entity_index]);
  }

  StateType get_alpha(const StateType& u) const
  {
    const auto alpha = implementation_->get_alpha(u)->first;
    StateType ret;
    std::copy(alpha.begin(), alpha.end(), ret.begin());
    return ret;
  }

  template <class FluxesMapType>
  void calculate_reconstructed_fluxes(const FieldVector<size_t, 3>& entity_indices,
                                      const FieldVector<bool, 3>& boundary_direction,
                                      FluxesMapType& precomputed_fluxes,
                                      const size_t dd) const
  {
    FieldVector<const QuadratureWeightsType*, 3> densities_stencil;
    for (size_t ii = 0; ii < 3; ++ii)
      if (entity_indices[ii] == size_t(-1))
        densities_stencil[ii] = &boundary_distribution_evaluations_[entity_indices[1]][dd][boundary_direction[ii]];
      else
        densities_stencil[ii] = &(*eta_ast_prime_evaluations_)[entity_indices[ii]];
    implementation_->template calculate_reconstructed_fluxes<slope, FluxesMapType>(
        densities_stencil, precomputed_fluxes, dd);
  }

  void apply_inverse_hessian(const size_t entity_index, const StateType& u, StateType& Hinv_u) const
  {
    implementation_->apply_inverse_hessian((*eta_ast_twoprime_evaluations_)[entity_index], u, Hinv_u);
  }

  void store_evaluations(const size_t entity_index, StateType& alpha, bool check = true)
  {
    implementation_->store_exp_evaluations(exp_evaluations_[entity_index], alpha);
    if (entropy != EntropyType::MaxwellBoltzmann) {
      implementation_->store_eta_ast_prime_vals(exp_evaluations_[entity_index], eta_ast_prime_storage_[entity_index]);
      implementation_->store_eta_ast_twoprime_vals(exp_evaluations_[entity_index],
                                                   eta_ast_twoprime_storage_[entity_index]);
    }
    set_eta_ast_pointers();
    // check for inf and nan and very low densities
    if (check) {
      const auto& u = get_u(entity_index);
      for (auto&& entry : u)
        if (std::isnan(entry) || std::isinf(entry))
          DUNE_THROW(Dune::MathError, "inf or nan in u!");
      const auto rho = basis_functions().density(u);
      if (rho < 1e-9) {
        alpha = basis_functions().alpha_iso(1e-8);
        store_evaluations(entity_index, alpha, false);
      }
    }
  }

  void set_eta_ast_pointers()
  {
    if (entropy == EntropyType::MaxwellBoltzmann) {
      eta_ast_prime_evaluations_ = &exp_evaluations_;
      eta_ast_twoprime_evaluations_ = &exp_evaluations_;
    } else {
      eta_ast_prime_evaluations_ = &eta_ast_prime_storage_;
      eta_ast_twoprime_evaluations_ = &eta_ast_twoprime_storage_;
    }
  }

  void store_boundary_evaluations(const std::function<RangeFieldType(const DomainType&)>& boundary_distribution,
                                  const size_t entity_index,
                                  const size_t intersection_index)
  {
    implementation_->store_boundary_distribution_evaluations(
        boundary_distribution_evaluations_[entity_index][intersection_index / 2][intersection_index % 2],
        boundary_distribution);
  }

  std::vector<QuadratureWeightsType>& exp_evaluations()
  {
    return exp_evaluations_;
  }

  const std::vector<QuadratureWeightsType>& exp_evaluations() const
  {
    return exp_evaluations_;
  }

  std::vector<QuadratureWeightsType>& eta_ast_prime_evaluations()
  {
    return *eta_ast_prime_evaluations_;
  }

  const std::vector<QuadratureWeightsType>& eta_ast_prime_evaluations() const
  {
    return *eta_ast_prime_evaluations_;
  }

  std::vector<QuadratureWeightsType>& eta_ast_twoprime_evaluations()
  {
    return *eta_ast_twoprime_evaluations_;
  }

  const std::vector<QuadratureWeightsType>& eta_ast_twoprime_evaluations() const
  {
    return *eta_ast_twoprime_evaluations_;
  }

  BoundaryQuadratureWeightsType& boundary_distribution_evaluations()
  {
    return boundary_distribution_evaluations_;
  }

  const BoundaryQuadratureWeightsType& boundary_distribution_evaluations() const
  {
    return boundary_distribution_evaluations_;
  }

private:
  std::shared_ptr<ImplementationType> implementation_;
  std::vector<QuadratureWeightsType> exp_evaluations_;
  std::vector<QuadratureWeightsType> eta_ast_prime_storage_;
  std::vector<QuadratureWeightsType> eta_ast_twoprime_storage_;
  std::vector<QuadratureWeightsType>* eta_ast_prime_evaluations_;
  std::vector<QuadratureWeightsType>* eta_ast_twoprime_evaluations_;
  BoundaryQuadratureWeightsType boundary_distribution_evaluations_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_KINETICCOORDS_HH
