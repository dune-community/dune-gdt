// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_LOCAL_INTEGRANDS_ESV2007_HH
#define DUNE_GDT_LOCAL_INTEGRANDS_ESV2007_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, to be used in the traits
template <class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType = void>
class LocalDiffusiveFluxEstimateESV2007Integrand;


namespace internal {


template <class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType = void>
class LocalDiffusiveFluxEstimateESV2007IntegrandTraits
{
  static_assert(XT::Functions::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be a localizable function.");
  static_assert(XT::Functions::is_localizable_function<DiffusiveFluxType>::value,
                "DiffusiveFluxType has to be a localizable function.");
  static_assert(XT::Functions::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be a localizable function.");
  static_assert(
      std::is_same<typename DiffusionFactorType::EntityType, typename DiffusiveFluxType::EntityType>::value
          && std::is_same<typename DiffusionFactorType::EntityType, typename DiffusionTensorType::EntityType>::value,
      "Types have to agree!");
  static_assert(
      std::is_same<typename DiffusionFactorType::DomainFieldType, typename DiffusiveFluxType::DomainFieldType>::value
          && std::is_same<typename DiffusionFactorType::DomainFieldType,
                          typename DiffusionTensorType::DomainFieldType>::value,
      "Types have to agree!");
  static_assert(DiffusionFactorType::dimDomain == DiffusiveFluxType::dimDomain
                    && DiffusionFactorType::dimDomain == DiffusionTensorType::dimDomain,
                "Dimensions have to agree");

public:
  typedef LocalDiffusiveFluxEstimateESV2007Integrand<DiffusionFactorType, DiffusiveFluxType, DiffusionTensorType>
      derived_type;
  typedef std::tuple<std::shared_ptr<typename DiffusionFactorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusionTensorType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusiveFluxType::LocalfunctionType>>
      LocalfunctionTupleType;
  typedef typename DiffusionFactorType::EntityType EntityType;
  typedef typename DiffusionFactorType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = DiffusionFactorType::dimDomain;
}; // class LocalDiffusiveFluxEstimateESV2007IntegrandTraits


template <class DiffusionType, class DiffusiveFluxType>
class LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionType, DiffusiveFluxType, void>
{
  static_assert(XT::Functions::is_localizable_function<DiffusionType>::value,
                "DiffusionType has to be a localizable function.");
  static_assert(XT::Functions::is_localizable_function<DiffusiveFluxType>::value,
                "DiffusiveFluxType has to be a localizable function.");
  static_assert(std::is_same<typename DiffusionType::EntityType, typename DiffusiveFluxType::EntityType>::value,
                "Types have to agree!");
  static_assert(
      std::is_same<typename DiffusionType::DomainFieldType, typename DiffusiveFluxType::DomainFieldType>::value,
      "Types have to agree!");
  static_assert(DiffusionType::dimDomain == DiffusiveFluxType::dimDomain, "Dimensions of domains have to agree");

public:
  typedef LocalDiffusiveFluxEstimateESV2007Integrand<DiffusionType, DiffusiveFluxType> derived_type;
  typedef std::tuple<std::shared_ptr<typename DiffusionType::LocalfunctionType>,
                     std::shared_ptr<typename DiffusiveFluxType::LocalfunctionType>>
      LocalfunctionTupleType;
  typedef typename DiffusionType::EntityType EntityType;
  typedef typename DiffusionType::DomainFieldType DomainFieldType;
  static const size_t dimDomain = DiffusionType::dimDomain;
}; // class LocalDiffusiveFluxEstimateESV2007IntegrandTraits


} // namespace internal


template <class DiffusionType, class DiffusiveFluxType>
class LocalDiffusiveFluxEstimateESV2007Integrand<DiffusionType, DiffusiveFluxType, void>
    : public LocalVolumeIntegrandInterface<internal::
                                               LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionType,
                                                                                                DiffusiveFluxType>,
                                           2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionType,
                                                                                                   DiffusiveFluxType>,
                                        2>
      BaseType;

public:
  typedef internal::LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionType, DiffusiveFluxType> Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  LocalDiffusiveFluxEstimateESV2007Integrand(const DiffusionType& diffusion, const DiffusiveFluxType& diffusive_flux)
    : diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
  {
  }

  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity), diffusive_flux_.local_function(entity));
  }

  /// \brief Extracts the local functions and calls the correct order method.
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType localFuncs,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion      = std::get<0>(localFuncs);
    const auto local_diffusive_flux = std::get<1>(localFuncs);
    return order(*local_diffusion, *local_diffusive_flux, testBase, ansatzBase);
  }

  /// \brief Extracts the local functions and calls the correct evaluate method.
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base,
                const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion      = std::get<0>(localFuncs);
    const auto local_diffusive_flux = std::get<1>(localFuncs);
    evaluate(*local_diffusion, *local_diffusive_flux, test_base, ansatz_base, local_point, ret);
  }

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rLD, size_t rCLD, size_t rLDF, size_t rCLDF, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLD, rCLD>& local_diffusion,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF>& local_diffusive_flux,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    // TODO: there is no way to guess the order of local_diffusion^-1, so we take local_diffusion.order()
    return local_diffusion.order()
           + (std::max(local_diffusion.order() + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()))
           + (std::max(local_diffusion.order() + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()));
  } // ... order(...)

  /// \}
  /// \name Actual implementation of evalaute.
  /// \{

  template <class R>
  void evaluate(
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain>& local_diffusive_flux,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& test_base,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& ansatz_base,
      const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    // evaluate local functions
    const auto diffusion_value       = local_diffusion.evaluate(local_point);
    const R one_over_diffusion_value = 1.0 / diffusion_value;
    const auto diffusive_flux_value  = local_diffusive_flux.evaluate(local_point);
    // evaluate test gradient
    const size_t rows         = test_base.size();
    const auto test_gradients = test_base.jacobian(local_point);
    auto left_sum             = test_gradients[0];
    assert(test_gradients.size() == rows);
    // evaluate ansatz gradient
    const size_t cols           = ansatz_base.size();
    const auto ansatz_gradients = ansatz_base.jacobian(local_point);
    auto right_sum              = ansatz_gradients[0];
    assert(ansatz_gradients.size() == rows);
    // compute
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      left_sum = test_gradients[ii];
      left_sum[0] *= diffusion_value;
      left_sum[0] += diffusive_flux_value;
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        right_sum = ansatz_gradients[jj];
        right_sum[0] *= diffusion_value;
        right_sum[0] += diffusive_flux_value;
        retRow[jj] = one_over_diffusion_value * (left_sum[0] * right_sum[0]);
      }
    }
  } // ... evaluate(...)

  /// \}

private:
  const DiffusionType& diffusion_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class LocalDiffusiveFluxEstimateESV2007Integrand< ..., void >


template <class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType>
class LocalDiffusiveFluxEstimateESV2007Integrand
    : public LocalVolumeIntegrandInterface<internal::
                                               LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionFactorType,
                                                                                                DiffusiveFluxType,
                                                                                                DiffusionTensorType>,
                                           2>
{
  typedef LocalVolumeIntegrandInterface<internal::LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionFactorType,
                                                                                                   DiffusiveFluxType,
                                                                                                   DiffusionTensorType>,
                                        2>
      BaseType;

public:
  typedef internal::LocalDiffusiveFluxEstimateESV2007IntegrandTraits<DiffusionFactorType, DiffusiveFluxType,
                                                                     DiffusionTensorType>
      Traits;
  using typename BaseType::LocalfunctionTupleType;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;

  LocalDiffusiveFluxEstimateESV2007Integrand(const DiffusionFactorType& diffusion_factor,
                                             const DiffusionTensorType& diffusion_tensor,
                                             const DiffusiveFluxType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {
  }

  /// \name Required by LocalVolumeIntegrandInterface< ..., 2 >
  /// \{

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  /// \brief Extracts the local functions and calls the correct order method.
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  size_t
  order(const LocalfunctionTupleType localFuncs,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& testBase,
        const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    const auto local_diffusive_flux   = std::get<2>(localFuncs);
    return order(*local_diffusion_factor, *local_diffusion_tensor, *local_diffusive_flux, testBase, ansatzBase);
  } // ... order(...)

  /// \brief Extracts the local functions and calls the correct evaluate method.
  template <class R, size_t rT, size_t rCT, size_t rA, size_t rCA>
  void evaluate(const LocalfunctionTupleType localFuncs,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
                const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base,
                const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    const auto local_diffusion_factor = std::get<0>(localFuncs);
    const auto local_diffusion_tensor = std::get<1>(localFuncs);
    const auto local_diffusive_flux   = std::get<2>(localFuncs);
    evaluate(*local_diffusion_factor,
             *local_diffusion_tensor,
             *local_diffusive_flux,
             test_base,
             ansatz_base,
             local_point,
             ret);
  } // ... evaluate(...)

  /// \}
  /// \name Actual implementation of order.
  /// \{

  template <class R, size_t rLD, size_t rCLD, size_t rLDT, size_t rCLDT, size_t rLDF, size_t rCLDF, size_t rT,
            size_t rCT, size_t rA, size_t rCA>
  size_t order(
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLD, rCLD>& local_diffusion_factor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDT, rCLDT>&
          local_diffusion_tensor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF>& local_diffusive_flux,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rT, rCT>& test_base,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, rA, rCA>& ansatz_base) const
  {
    // TODO: there is no way to guess the order of (local_diffusion_factor * local_diffusion_tensor)^-1,
    //       so we take local_diffusion_factor.order() + local_diffusion_tensor.order()
    const size_t local_diffusion_order = local_diffusion_factor.order() + local_diffusion_tensor.order();
    return local_diffusion_order
           + (std::max(local_diffusion_order + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()))
           + (std::max(local_diffusion_order + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                       local_diffusive_flux.order()));
  } // ... order(...)

  /// \}
  /// \name Actual implementation of evalaute.
  /// \{

  template <class R>
  void evaluate(
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, 1>& local_diffusion_factor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain>&
          local_diffusion_tensor,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, R, dimDomain>& local_diffusive_flux,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& test_base,
      const XT::Functions::LocalfunctionSetInterface<EntityType, DomainFieldType, dimDomain, R, 1>& ansatz_base,
      const Dune::FieldVector<DomainFieldType, dimDomain>& local_point, Dune::DynamicMatrix<R>& ret) const
  {
    typedef FieldVector<R, dimDomain> DomainType;
    DomainType left_sum(0);
    DomainType right_sum(0);
    // evaluate local functions
    const auto diffusion_factor_value = local_diffusion_factor.evaluate(local_point);
    typedef XT::Common::FieldMatrix<R, dimDomain, dimDomain> TensorType;
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point);
    const TensorType diffusion_value        = diffusion_tensor_value * diffusion_factor_value;
    // TODO: there is no documented way to assert that the inversion was successfull, so find one or check the matrix
    //       beforehand
    TensorType one_over_diffusion_value = diffusion_tensor_value;
    one_over_diffusion_value.invert();
    one_over_diffusion_value /= diffusion_factor_value;
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
        retRow[jj] = (one_over_diffusion_value * left_sum) * right_sum;
      }
    }
  } // ... evaluate(...)

  /// \}

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class LocalDiffusiveFluxEstimateESV2007Integrand


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_INTEGRANDS_ESV2007_HH
