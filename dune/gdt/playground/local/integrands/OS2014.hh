// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_PLAYGROUND_LOCAL_INTEGRANDS_OS2014_HH
#define DUNE_GDT_PLAYGROUND_LOCAL_INTEGRANDS_OS2014_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/local/integrands/interfaces.hh>

namespace Dune {
namespace GDT {


// forward, to be used in the traits
template <class DiffusionFactorImp, class DiffusionFactorHatImp, class DiffusionTensorImp, class DiffusiveFluxImp>
class LocalDiffusiveFluxEstimateStarOS2014Integrand;


namespace internal {


template <class DiffusionFactorType, class DiffusionFactorHatType, class DiffusionTensorType, class DiffusiveFluxType>
class LocalDiffusiveFluxEstimateStarOS2014IntegrandTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be a localizable function.");
  static_assert(Stuff::is_localizable_function<DiffusionFactorHatType>::value,
                "DiffusionFactorHatType has to be a localizable function.");
  static_assert(Stuff::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be a localizable function.");
  static_assert(Stuff::is_localizable_function<DiffusiveFluxType>::value,
                "DiffusiveFluxType has to be a localizable function.");
  static_assert(
      std::is_same<typename DiffusionFactorType::EntityType, typename DiffusionTensorType::EntityType>::value
          && std::is_same<typename DiffusionFactorHatType::EntityType, typename DiffusiveFluxType::EntityType>::value
          && std::is_same<typename DiffusionFactorType::EntityType, typename DiffusionFactorHatType::EntityType>::value,
      "Types have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorType::DomainFieldType, typename DiffusionTensorType::DomainFieldType>::value
          && std::is_same<typename DiffusionFactorHatType::DomainFieldType,
                          typename DiffusiveFluxType::DomainFieldType>::value
          && std::is_same<typename DiffusionFactorType::DomainFieldType,
                          typename DiffusionFactorHatType::DomainFieldType>::value,
      "Types have to agree!");
  static_assert(DiffusionFactorType::dimDomain == DiffusionTensorType::dimDomain
                    && DiffusionFactorHatType::dimDomain == DiffusiveFluxType::dimDomain
                    && DiffusionFactorType::dimDomain == DiffusionFactorHatType::dimDomain,
                "Dimensions have to agree");

public:
  typedef LocalDiffusiveFluxEstimateStarOS2014Integrand<DiffusionFactorType, DiffusionFactorHatType,
                                                        DiffusionTensorType, DiffusiveFluxType>
      derived_type;
  typedef std::tuple<std::shared_ptr<typename DiffusionFactorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionFactorHatType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionTensorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusiveFluxType::LocalfunctionType>>
      LocalfunctionTupleType;
  typedef typename DiffusionFactorType::EntityType EntityType;
  typedef typename DiffusionFactorType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = DiffusionFactorType::dimDomain;
}; // class DiffusiveFluxEstimateTraits


} // namespace internal


template <class DiffusionFactorType, class DiffusionFactorHatType, class DiffusionTensorType, class DiffusiveFluxType>
class LocalDiffusiveFluxEstimateStarOS2014Integrand
    : public LocalVolumeIntegrandInterface<internal::
                                               LocalDiffusiveFluxEstimateStarOS2014IntegrandTraits<DiffusionFactorType,
                                                                                                   DiffusionFactorHatType,
                                                                                                   DiffusionTensorType,
                                                                                                   DiffusiveFluxType>,
                                           2>
{
  typedef LocalVolumeIntegrandInterface<internal::
                                            LocalDiffusiveFluxEstimateStarOS2014IntegrandTraits<DiffusionFactorType,
                                                                                                DiffusionFactorHatType,
                                                                                                DiffusionTensorType,
                                                                                                DiffusiveFluxType>,
                                        2>
      BaseType;

public:
  typedef internal::LocalDiffusiveFluxEstimateStarOS2014IntegrandTraits<DiffusionFactorType, DiffusionFactorHatType,
                                                                        DiffusionTensorType, DiffusiveFluxType>
      Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  LocalDiffusiveFluxEstimateStarOS2014Integrand(const DiffusionFactorType& diffusion_factor,
                                                const DiffusionFactorHatType& diffusion_factor_hat,
                                                const DiffusionTensorType& diffusion_tensor,
                                                const DiffusiveFluxType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_factor_hat_(diffusion_factor_hat)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {
  }

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_factor_hat_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType& localFuncs,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion_factor     = std::get<0>(localFuncs);
    const auto local_diffusion_factor_hat = std::get<1>(localFuncs);
    const auto local_diffusion_tensor     = std::get<2>(localFuncs);
    const auto local_diffusive_flux       = std::get<3>(localFuncs);
    return order(*local_diffusion_factor,
                 *local_diffusion_factor_hat,
                 *local_diffusion_tensor,
                 *local_diffusive_flux,
                 testBase,
                 ansatzBase);
  } // ... order(...)

  template <class R, size_t rLD, size_t rCLD, size_t rLDT, size_t rCLDT, size_t rLDF, size_t rCLDF, size_t rT,
            size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLD, rCLD>& local_diffusion_factor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLD, rCLD>&
          local_diffusion_factor_hat,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDT, rCLDT>&
          local_diffusion_tensor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF>& local_diffusive_flux,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    // TODO: there is no way to guess the order of (local_diffusion_factor * local_diffusion_tensor)^-1,
    //       so we take local_diffusion_factor.order() + local_diffusion_tensor.order()
    const size_t local_diffusion_order     = local_diffusion_factor.order() + local_diffusion_tensor.order();
    const size_t local_diffusion_hat_order = local_diffusion_factor_hat.order() + local_diffusion_tensor.order();
    return local_diffusion_hat_order
           + (std::max(local_diffusion_order + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()))
           + (std::max(local_diffusion_order + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()));
  } // ... order(...)

  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
                const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base,
                const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor     = std::get<0>(localFuncs);
    const auto local_diffusion_factor_hat = std::get<1>(localFuncs);
    const auto local_diffusion_tensor     = std::get<2>(localFuncs);
    const auto local_diffusive_flux       = std::get<3>(localFuncs);
    evaluate(*local_diffusion_factor,
             *local_diffusion_factor_hat,
             *local_diffusion_tensor,
             *local_diffusive_flux,
             test_base,
             ansatz_base,
             local_point,
             ret);
  } // ... evaluate(...)

  template <class R>
  void evaluate(
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor_hat,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          local_diffusion_tensor,
      const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain>& local_diffusive_flux,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& test_base,
      const Stuff::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& ansatz_base,
      const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    typedef FieldVector<R, dimDomain> DomainType;
    DomainType left_sum(0);
    DomainType right_sum(0);
    // evaluate local functions
    const auto diffusion_factor_value     = local_diffusion_factor.evaluate(local_point);
    const auto diffusion_factor_hat_value = local_diffusion_factor_hat.evaluate(local_point);
    typedef Stuff::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point);
    const TensorType diffusion_value        = diffusion_tensor_value * diffusion_factor_value;
    // TODO: there is no documented way to assert that the inversion was successfull, so find one or check the matrix
    //       beforehand
    TensorType one_over_diffusion_hat_value = diffusion_tensor_value;
    one_over_diffusion_hat_value.invert();
    one_over_diffusion_hat_value /= diffusion_factor_hat_value;
    const auto diffusive_flux_value = local_diffusive_flux.evaluate(local_point);
    // evaluate test gradient
    const size_t rows         = test_base.size();
    const auto test_gradients = test_base.jacobian(local_point);
    assert(test_gradients.size() == rows);
    // evaluate ansatz gradient
    const size_t cols           = ansatz_base.size();
    const auto ansatz_gradients = ansatz_base.jacobian(local_point);
    assert(ansatz_gradients.size() == rows);
    // compute
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      left_sum = diffusion_value * test_gradients[ii];
      left_sum += diffusive_flux_value;
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        right_sum = diffusion_value * ansatz_gradients[jj];
        right_sum += diffusive_flux_value;
        retRow[jj] = (one_over_diffusion_hat_value * left_sum) * right_sum;
      }
    }
  } // ... evaluate(...)

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionFactorHatType& diffusion_factor_hat_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class LocalDiffusiveFluxEstimateStarOS2014Integrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_LOCAL_INTEGRANDS_OS2014_HH
