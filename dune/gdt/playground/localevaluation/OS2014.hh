// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALEVALUATION_OS2014_HH
#define DUNE_GDT_LOCALEVALUATION_OS2014_HH

#include <tuple>
#include <memory>

#include <dune/stuff/common/disable_warnings.hh>
# include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "../../localevaluation/interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace OS2014 {


// forward, to be used in the traits
template< class DiffusionFactorType, class DiffusionFactorHatType, class DiffusionTensorType, class DiffusiveFluxType >
class DiffusiveFluxEstimateStar;


namespace internal {


template< class DiffusionFactorType, class DiffusionFactorHatType, class DiffusionTensorType, class DiffusiveFluxType >
class DiffusiveFluxEstimateStarTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorType >::value,
                "DiffusionFactorType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorHatType >::value,
                "DiffusionFactorHatType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionTensorType >::value,
                "DiffusionTensorType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxType >::value,
                "DiffusiveFluxType has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef DiffusiveFluxEstimateStar
      < DiffusionFactorType, DiffusionFactorHatType, DiffusionTensorType, DiffusiveFluxType > derived_type;
}; // class DiffusiveFluxEstimateTraits


} // namespace internal


template< class DiffusionFactorType, class DiffusionFactorHatType, class DiffusionTensorType, class DiffusiveFluxType >
class DiffusiveFluxEstimateStar
  : public LocalEvaluation::Codim0Interface< internal::DiffusiveFluxEstimateStarTraits< DiffusionFactorType,
                                                                                        DiffusionFactorHatType,
                                                                                        DiffusionTensorType,
                                                                                        DiffusiveFluxType >,
                                             2 >
{
public:
  typedef internal::DiffusiveFluxEstimateStarTraits
      < DiffusionFactorType, DiffusionFactorHatType, DiffusionTensorType, DiffusiveFluxType > Traits;

  DiffusiveFluxEstimateStar(const DiffusionFactorType& diffusion_factor,
                            const DiffusionFactorHatType& diffusion_factor_hat,
                            const DiffusionTensorType& diffusion_tensor,
                            const DiffusiveFluxType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_factor_hat_(diffusion_factor_hat)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType    LocalDiffusionFactorType;
    typedef typename DiffusionFactorHatType::LocalfunctionType LocalDiffusionFactorHatType;
    typedef typename DiffusionTensorType::LocalfunctionType    LocalDiffusionTensorType;
    typedef typename DiffusiveFluxType::LocalfunctionType      LocalDiffusiveFluxType;
  public:
    typedef std::tuple< std::shared_ptr< LocalDiffusionFactorType >,
                        std::shared_ptr< LocalDiffusionFactorHatType >,
                        std::shared_ptr< LocalDiffusionTensorType >,
                        std::shared_ptr< LocalDiffusiveFluxType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_factor_hat_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase)
  {
    const auto local_diffusion_factor     = std::get< 0 >(localFuncs);
    const auto local_diffusion_factor_hat = std::get< 1 >(localFuncs);
    const auto local_diffusion_tensor     = std::get< 2 >(localFuncs);
    const auto local_diffusive_flux       = std::get< 3 >(localFuncs);
    return redirect_order(*local_diffusion_factor,
                          *local_diffusion_factor_hat,
                          *local_diffusion_tensor,
                          *local_diffusive_flux,
                          testBase, ansatzBase);
  } // ... order(...)

private:
  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static size_t redirect_order(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& local_diffusion_factor,
                               const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& local_diffusion_factor_hat,
                               const Stuff::LocalfunctionInterface< E, D, d, R, rLDT, rCLDT >& local_diffusion_tensor,
                               const Stuff::LocalfunctionInterface< E, D, d, R, rLDF, rCLDF >& local_diffusive_flux,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base)
  {
    // TODO: there is no way to guess the order of (local_diffusion_factor * local_diffusion_tensor)^-1,
    //       so we take local_diffusion_factor.order() + local_diffusion_tensor.order()
    const size_t local_diffusion_order     = local_diffusion_factor.order()     + local_diffusion_tensor.order();
    const size_t local_diffusion_hat_order = local_diffusion_factor_hat.order() + local_diffusion_tensor.order();
    return local_diffusion_hat_order
        + (std::max(local_diffusion_order + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()))
        + (std::max(local_diffusion_order + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()));
  } // ... redirect_order(...)

public:
  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base,
                       const Dune::FieldVector< D, d >& local_point,
                       Dune::DynamicMatrix< R >& ret)
  {
    const auto local_diffusion_factor     = std::get< 0 >(localFuncs);
    const auto local_diffusion_factor_hat = std::get< 1 >(localFuncs);
    const auto local_diffusion_tensor     = std::get< 2 >(localFuncs);
    const auto local_diffusive_flux       = std::get< 3 >(localFuncs);
    redirect_evaluate(*local_diffusion_factor,
                      *local_diffusion_factor_hat,
                      *local_diffusion_tensor,
                      *local_diffusive_flux,
                      test_base, ansatz_base, local_point, ret);
  } // ... evaluate(...)

private:
  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDh, int rCLDh, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static void redirect_evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& /*local_diffusion_factor*/,
                                const Stuff::LocalfunctionInterface< E, D, d, R, rLDh, rCLDh >& /*local_diffusion_factor_hat*/,
                                const Stuff::LocalfunctionInterface< E, D, d, R, rLDT, rCLDT >& /*local_diffusion_tensor*/,
                                const Stuff::LocalfunctionInterface< E, D, d, R, rLDF, rCLDF >& /*local_diffusive_flux*/,
                                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*test_base*/,
                                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatz_base*/,
                                const Dune::FieldVector< D, d >& /*local_point*/,
                                Dune::DynamicMatrix< R >& /*ret*/)
  {
    static_assert(AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  template< class E, class D, int d, class R >
  static void redirect_evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1 >& local_diffusion_factor,
                                const Stuff::LocalfunctionInterface< E, D, d, R, 1 >& local_diffusion_factor_hat,
                                const Stuff::LocalfunctionInterface< E, D, d, R, d, d >& local_diffusion_tensor,
                                const Stuff::LocalfunctionInterface< E, D, d, R, d >& local_diffusive_flux,
                                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1 >& test_base,
                                const Stuff::LocalfunctionSetInterface< E, D, d, R, 1 >& ansatz_base,
                                const Dune::FieldVector< D, d >& local_point,
                                Dune::DynamicMatrix< R >& ret)
  {
    typedef FieldVector< R, d > DomainType;
    DomainType left_sum(0);
    DomainType right_sum(0);
    // evaluate local functions
    const auto diffusion_factor_value     = local_diffusion_factor.evaluate(local_point);
    const auto diffusion_factor_hat_value = local_diffusion_factor_hat.evaluate(local_point);
    typedef Stuff::Common::FieldMatrix< R, d, d > TensorType;
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point);
    const TensorType diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // TODO: there is no documented way to assert that the inversion was successfull, so find one or check the matrix
    //       beforehand
    TensorType one_over_diffusion_hat_value = diffusion_tensor_value;
    one_over_diffusion_hat_value.invert();
    one_over_diffusion_hat_value /= diffusion_factor_hat_value;
    const auto diffusive_flux_value = local_diffusive_flux.evaluate(local_point);
    // evaluate test gradient
    const size_t rows = test_base.size();
    const auto test_gradients = test_base.jacobian(local_point);
    assert(test_gradients.size() == rows);
    // evaluate ansatz gradient
    const size_t cols = ansatz_base.size();
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
  } // ... redirect_evaluate(...)

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionFactorHatType& diffusion_factor_hat_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class DiffusiveFluxEstimateStar


} // namespace OS2014
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_OS2014_HH
