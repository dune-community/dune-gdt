// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2016 - 2017)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
#define DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH

#include <list>
#include <memory>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/vector_less.hh>

#include <dune/xt/functions/interfaces/flux-function.hh>

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux_implementations.hh>

namespace Dune {
namespace GDT {


// Caches a specified number of (u, alpha) pairs. If the cache is full and another pair is added, the oldest existing
// pair is dropped.
template <class KeyVectorType, class ValueVectorType>
class EntropyLocalCache
{
public:
  using MapType = typename std::map<KeyVectorType, ValueVectorType, XT::Common::VectorLess>;
  using IteratorType = typename MapType::iterator;
  using ConstIteratorType = typename MapType::const_iterator;
  using RangeFieldType = typename XT::Common::VectorAbstraction<KeyVectorType>::ScalarType;

  EntropyLocalCache(const size_t capacity = 0)
    : capacity_(capacity)
  {}

  void insert(const KeyVectorType& u, const ValueVectorType& alpha)
  {
    cache_.insert(std::make_pair(u, alpha));
    keys_.push_back(u);
    if (cache_.size() > capacity_) {
      cache_.erase(keys_.front());
      keys_.pop_front();
    }
  }

  std::pair<RangeFieldType, ConstIteratorType> find_closest(const KeyVectorType& u) const
  {
    ConstIteratorType ret = cache_.begin();
    if (ret == end())
      return std::make_pair(std::numeric_limits<RangeFieldType>::max(), ret);
    auto diff = u - ret->first;
    // use infinity_norm as distance
    RangeFieldType distance = infinity_norm(diff);
    auto it = ret;
    while (++it != end()) {
      if (XT::Common::FloatCmp::eq(distance, 0.))
        break;
      diff = u - it->first;
      RangeFieldType new_distance = infinity_norm(diff);
      if (new_distance < distance) {
        distance = new_distance;
        ret = it;
      }
    }
    return std::make_pair(distance, ret);
  }

  IteratorType begin()
  {
    return cache_.begin();
  }

  ConstIteratorType begin() const
  {
    return cache_.begin();
  }

  IteratorType end()
  {
    return cache_.end();
  }

  ConstIteratorType end() const
  {
    return cache_.end();
  }

private:
  static RangeFieldType infinity_norm(const KeyVectorType& vec)
  {
    RangeFieldType ret = std::abs(vec[0]);
    for (size_t ii = 1; ii < vec.size(); ++ii)
      ret = std::max(ret, std::abs(vec[ii]));
    return ret;
  }

  size_t capacity_;
  MapType cache_;
  std::list<KeyVectorType> keys_;
};


// This flux function only does the caching, which is used for providing a good initial guess and to avoid solving the
// same optimization problem twice. The actual implementations of the optimization algorithm are in
// entropyflux_implementations.hh
template <class GridViewImp, class MomentBasisImp>
class EntropyBasedFluxFunction
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
  using ThisType = EntropyBasedFluxFunction;

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
  using LocalCacheType = EntropyLocalCache<StateType, VectorType>;
  static const size_t cache_size = 4 * dimFlux + 2;

  explicit EntropyBasedFluxFunction(
      const GridViewType& grid_view,
      const MomentBasis& basis_functions,
      const RangeFieldType tau = 1e-9,
      const bool disable_realizability_check = false,
      const RangeFieldType epsilon_gamma = 0.01,
      const RangeFieldType chi = 0.5,
      const RangeFieldType xi = 1e-3,
      const std::vector<RangeFieldType> r_sequence = {0, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 5e-2, 0.1, 0.5, 1},
      const size_t k_0 = 500,
      const size_t k_max = 1000,
      const RangeFieldType epsilon = std::pow(2, -52))
    : index_set_(grid_view.indexSet())
    , use_thread_cache_(true)
    , use_entity_cache_(true)
    , entity_caches_(index_set_.size(0), LocalCacheType(cache_size))
    , mutexes_(index_set_.size(0))
    , implementation_(std::make_shared<ImplementationType>(
          basis_functions, tau, disable_realizability_check, epsilon_gamma, chi, xi, r_sequence, k_0, k_max, epsilon))
  {}

  void enable_thread_cache()
  {
    use_thread_cache_ = true;
  }

  void disable_thread_cache()
  {
    use_thread_cache_ = false;
  }

  void enable_entity_cache()
  {
    use_entity_cache_ = true;
  }

  void disable_entity_cache()
  {
    use_entity_cache_ = false;
  }

  static const constexpr bool available = true;

  class Localfunction : public LocalFunctionType
  {
    using BaseType = LocalFunctionType;

  public:
    using typename BaseType::DynamicJacobianRangeType;
    using typename BaseType::E;
    using typename BaseType::RangeReturnType;

    Localfunction(const IndexSetType& index_set,
                  std::vector<LocalCacheType>& entity_caches,
                  const bool use_thread_cache,
                  const bool use_entity_cache,
                  std::vector<std::mutex>& mutexes,
                  const ImplementationType& implementation)
      : index_set_(index_set)
      , thread_cache_(cache_size)
      , entity_caches_(entity_caches)
      , use_thread_cache_(use_thread_cache)
      , use_entity_cache_(use_entity_cache)
      , mutexes_(mutexes)
      , implementation_(implementation)
    {}

    void post_bind(const E& element) override final
    {
      const auto index = index_set_.index(element);
      entity_cache_ = &(entity_caches_[index]);
      mutex_ = &(mutexes_[index]);
    }

    int order(const XT::Common::Parameter&) const override final
    {
      return 1.;
    }

    std::unique_ptr<AlphaReturnType> get_alpha(const StateType& u, const bool regularize) const
    {
      // find starting point. Candidates: alpha_iso and the entries in the two caches
      std::lock_guard<std::mutex> DUNE_UNUSED(guard)(*mutex_);
      const auto& basis_functions = implementation_.basis_functions();
      static const auto u_iso = basis_functions.u_iso();
      const auto density = basis_functions.density(u);
      const auto alpha_iso = basis_functions.alpha_iso(density);
      const auto u_iso_scaled = u_iso * density;
      // calculate (inf-norm) distance to isotropic moment with same density
      RangeFieldType distance = (u - u_iso_scaled).infinity_norm();
      VectorType alpha_start = XT::Common::convert_to<VectorType>(alpha_iso);
      if (!XT::Common::FloatCmp::eq(distance, 0.) && use_entity_cache_) {
        // calculate distance to closest moment in entity_cache
        const auto entity_cache_dist_and_it = entity_cache_->find_closest(u);
        const auto& entity_cache_dist = entity_cache_dist_and_it.first;
        if (entity_cache_dist < distance) {
          distance = entity_cache_dist;
          alpha_start = entity_cache_dist_and_it.second->second;
        }
        if (!XT::Common::FloatCmp::eq(distance, 0.) && use_thread_cache_) {
          // calculate distance to closest moment in thread_cache
          const auto thread_cache_dist_and_it = thread_cache_.find_closest(u);
          const auto& thread_cache_dist = thread_cache_dist_and_it.first;
          if (thread_cache_dist < distance) {
            distance = thread_cache_dist;
            alpha_start = thread_cache_dist_and_it.second->second;
          }
        }
      }
      // If alpha_start is already the solution, we are finished. Else start optimization.
      if (XT::Common::FloatCmp::eq(distance, 0.)) {
        return std::make_unique<AlphaReturnType>(std::make_pair(alpha_start, std::make_pair(u, 0.)));
      } else {
        auto ret = implementation_.get_alpha(u, alpha_start, regularize);
        if (use_entity_cache_)
          entity_cache_->insert(ret->second.first, ret->first);
        if (use_thread_cache_)
          thread_cache_.insert(ret->second.first, ret->first);
        return ret;
      }
    }

    virtual RangeReturnType evaluate(const DomainType& /*point_in_reference_element*/,
                                     const StateType& u,
                                     const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      const auto alpha = get_alpha(u, true)->first;
      return implementation_.evaluate_with_alpha(alpha);
    }

    virtual void jacobian(const DomainType& /*point_in_reference_element*/,
                          const StateType& u,
                          DynamicJacobianRangeType& result,
                          const XT::Common::Parameter& /*param*/ = {}) const override final
    {
      const auto alpha = get_alpha(u, true)->first;
      implementation_.jacobian_with_alpha(alpha, result);
    } // ... jacobian(...)

  private:
    const IndexSetType& index_set_;
    mutable LocalCacheType thread_cache_;
    std::vector<LocalCacheType>& entity_caches_;
    const bool use_thread_cache_;
    const bool use_entity_cache_;
    std::vector<std::mutex>& mutexes_;
    const ImplementationType& implementation_;
    LocalCacheType* entity_cache_;
    std::mutex* mutex_;
  }; // class Localfunction

  bool x_dependent() const override final
  {
    return false;
  }

  std::unique_ptr<LocalFunctionType> local_function() const override final
  {
    return std::make_unique<Localfunction>(
        index_set_, entity_caches_, use_thread_cache_, use_entity_cache_, mutexes_, *implementation_);
  }

  virtual std::unique_ptr<Localfunction> derived_local_function() const
  {
    return std::make_unique<Localfunction>(
        index_set_, entity_caches_, use_thread_cache_, use_entity_cache_, mutexes_, *implementation_);
  }

  template <class StateTp, class RetType>
  void evaluate_kinetic_flux(const E& inside_entity,
                             const E& outside_entity,
                             const StateTp& u_i,
                             const StateTp& u_j,
                             const DomainType& n_ij,
                             const size_t dd,
                             RetType& ret) const
  {
    // calculate \sum_{i=1}^d < \omega_i m G_\alpha(u) > n_i
    const auto local_func = derived_local_function();
    local_func->bind(inside_entity);
    const auto alpha_i = local_func->get_alpha(u_i, true)->first;
    local_func->bind(outside_entity);
    const auto alpha_j = local_func->get_alpha(u_j, true)->first;
    ret = implementation_->evaluate_kinetic_flux_with_alphas(alpha_i, alpha_j, n_ij, dd);
  } // StateType evaluate_kinetic_flux(...)


  // Returns alpha(u), starting from alpha_iso. To get better performance when calculating several alphas, use
  // Localfunction's get_alpha
  std::unique_ptr<AlphaReturnType> get_alpha(const StateType& u, const bool regularize) const
  {
    const auto& basis_functions = implementation_->basis_functions();
    const auto density = basis_functions.density(u);
    auto alpha_iso = implementation_->get_isotropic_alpha(density);
    return implementation_->get_alpha(u, *alpha_iso, regularize);
  }

  const MomentBasis& basis_functions() const
  {
    return implementation_->basis_functions();
  }

  const IndexSetType& index_set_;
  bool use_thread_cache_;
  bool use_entity_cache_;
  mutable std::vector<LocalCacheType> entity_caches_;
  mutable std::vector<std::mutex> mutexes_;
  std::shared_ptr<ImplementationType> implementation_;
};

template <class GridViewImp, class MomentBasisImp>
const size_t EntropyBasedFluxFunction<GridViewImp, MomentBasisImp>::cache_size;


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_ENTROPYBASED_HH
