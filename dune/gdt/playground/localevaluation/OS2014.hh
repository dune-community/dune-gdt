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
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/functions/interfaces.hh>

#include "../../localevaluation/interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace OS2014 {


// forward, to be used in the traits
template< class DiffusionFactorImp, class DiffusionFactorHatImp, class DiffusionTensorImp, class DiffusiveFluxImp >
class DiffusiveFluxEstimateStar;


namespace internal {


template< class DiffusionFactorImp, class DiffusionFactorHatImp, class DiffusionTensorImp, class DiffusiveFluxImp >
class DiffusiveFluxEstimateStarTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorImp >::value,
                "DiffusionFactorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorHatImp >::value,
                "DiffusionFactorHatImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionTensorImp >::value,
                "DiffusionTensorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxImp >::value,
                "DiffusiveFluxImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_same< typename DiffusionFactorImp::EntityType,
                              typename DiffusionTensorImp::EntityType >::value &&
                std::is_same< typename DiffusionFactorHatImp::EntityType,
                              typename DiffusiveFluxImp::EntityType >::value &&
                std::is_same< typename DiffusionFactorImp::EntityType,
                              typename DiffusionFactorHatImp::EntityType >::value,
                "EntityTypes have to agree!");
  static_assert(std::is_same< typename DiffusionFactorImp::DomainFieldType,
                              typename DiffusionTensorImp::DomainFieldType >::value &&
                std::is_same< typename DiffusionFactorHatImp::DomainFieldType,
                              typename DiffusiveFluxImp::DomainFieldType >::value &&
                std::is_same< typename DiffusionFactorImp::DomainFieldType,
                              typename DiffusionFactorHatImp::DomainFieldType >::value,
                "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain &&
                DiffusionFactorHatImp::dimDomain == DiffusiveFluxImp::dimDomain &&
                DiffusionFactorImp::dimDomain == DiffusionFactorHatImp::dimDomain,
                "Dimensions of domains have to agree");

public:
  typedef DiffusiveFluxEstimateStar
      < DiffusionFactorImp, DiffusionFactorHatImp, DiffusionTensorImp, DiffusiveFluxImp > derived_type;

  typedef DiffusionFactorImp                                                  LocalizableDiffusionFactorFunctionType;
  typedef DiffusionFactorHatImp                                               LocalizableDiffusionFactorHatFunctionType;
  typedef DiffusionTensorImp                                                  LocalizableDiffusionTensorFunctionType;
  typedef DiffusiveFluxImp                                                    LocalizableDiffusiveFluxFunctionType;
  typedef typename LocalizableDiffusionFactorFunctionType::LocalfunctionType      LocalDiffusionFactorFunctionType;
  typedef typename LocalizableDiffusionFactorHatFunctionType::LocalfunctionType   LocalDiffusionFactorHatFunctionType;
  typedef typename LocalizableDiffusionTensorFunctionType::LocalfunctionType      LocalDiffusionTensorFunctionType;
  typedef typename LocalizableDiffusiveFluxFunctionType::LocalfunctionType        LocalDiffusiveFluxFunctionType;

  typedef std::tuple< std::shared_ptr< LocalDiffusionFactorFunctionType >,
                      std::shared_ptr< LocalDiffusionFactorHatFunctionType >,
                      std::shared_ptr< LocalDiffusionTensorFunctionType >,
                      std::shared_ptr< LocalDiffusiveFluxFunctionType > >         LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorFunctionType::EntityType             EntityType;
  typedef typename LocalizableDiffusionFactorFunctionType::DomainFieldType        DomainFieldType;
  static const unsigned int dimDomain = LocalizableDiffusionFactorFunctionType::dimDomain;
}; // class DiffusiveFluxEstimateTraits


} // namespace internal


template< class DiffusionFactorImp, class DiffusionFactorHatImp, class DiffusionTensorImp, class DiffusiveFluxImp >
class DiffusiveFluxEstimateStar
  : public LocalEvaluation::Codim0Interface< internal::DiffusiveFluxEstimateStarTraits< DiffusionFactorImp,
                                                                                        DiffusionFactorHatImp,
                                                                                        DiffusionTensorImp,
                                                                                        DiffusiveFluxImp >,
                                             2 >
{
public:
  typedef internal::DiffusiveFluxEstimateStarTraits
      < DiffusionFactorImp, DiffusionFactorHatImp, DiffusionTensorImp, DiffusiveFluxImp > Traits;
  typedef typename Traits::LocalizableDiffusionFactorFunctionType             LocalizableDiffusionFactorFunctionType;
  typedef typename Traits::LocalizableDiffusionFactorHatFunctionType          LocalizableDiffusionFactorHatFunctionType;
  typedef typename Traits::LocalizableDiffusionTensorFunctionType             LocalizableDiffusionTensorFunctionType;
  typedef typename Traits::LocalizableDiffusiveFluxFunctionType               LocalizableDiffusiveFluxFunctionType;
  typedef typename Traits::LocalfunctionTupleType                             LocalfunctionTupleType;
  typedef typename Traits::EntityType                                         EntityType;
  typedef typename Traits::DomainFieldType                                    DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  DiffusiveFluxEstimateStar(const LocalizableDiffusionFactorFunctionType& diffusion_factor,
                            const LocalizableDiffusionFactorHatFunctionType& diffusion_factor_hat,
                            const LocalizableDiffusionTensorFunctionType& diffusion_tensor,
                            const LocalizableDiffusiveFluxFunctionType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_factor_hat_(diffusion_factor_hat)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_factor_hat_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalfunctionTupleType& localFuncs,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase) const
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
  template< class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  size_t redirect_order(const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& local_diffusion_factor,
                        const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& local_diffusion_factor_hat,
                        const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLDT, rCLDT >& local_diffusion_tensor,
                        const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF >& local_diffusive_flux,
                        const Stuff::LocalfunctionSetInterface
                            < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                        const Stuff::LocalfunctionSetInterface
                            < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base) const
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
  template< class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalfunctionTupleType& localFuncs,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
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
  template< class R, int rLD, int rCLD, int rLDh, int rCLDh, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& /*local_diffusion_factor*/,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLDh, rCLDh >& /*local_diffusion_factor_hat*/,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLDT, rCLDT >& /*local_diffusion_tensor*/,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF >& /*local_diffusive_flux*/,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*test_base*/,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatz_base*/,
                         const Dune::FieldVector< DomainFieldType, dimDomain >& /*local_point*/,
                         Dune::DynamicMatrix< R >& /*ret*/) const
  {
    static_assert(AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  template< class R >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& local_diffusion_factor,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& local_diffusion_factor_hat,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, dimDomain, dimDomain >& local_diffusion_tensor,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, dimDomain >& local_diffusive_flux,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& test_base,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& ansatz_base,
                         const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                         Dune::DynamicMatrix< R >& ret) const
  {
    typedef FieldVector< R, dimDomain > DomainType;
    DomainType left_sum(0);
    DomainType right_sum(0);
    // evaluate local functions
    const auto diffusion_factor_value     = local_diffusion_factor.evaluate(local_point);
    const auto diffusion_factor_hat_value = local_diffusion_factor_hat.evaluate(local_point);
    typedef Stuff::Common::FieldMatrix< R, dimDomain, dimDomain > TensorType;
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

  const LocalizableDiffusionFactorFunctionType& diffusion_factor_;
  const LocalizableDiffusionFactorHatFunctionType& diffusion_factor_hat_;
  const LocalizableDiffusionTensorFunctionType& diffusion_tensor_;
  const LocalizableDiffusiveFluxFunctionType& diffusive_flux_;
}; // class DiffusiveFluxEstimateStar


} // namespace OS2014
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_OS2014_HH
