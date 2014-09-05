// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALEVALUATION_ESV2007_HH
#define DUNE_GDT_LOCALEVALUATION_ESV2007_HH

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
namespace ESV2007 {


// forward, to be used in the traits
template< class DiffusionFactorImp, class DiffusiveFluxImp, class DiffusionTensorImp = void >
class DiffusiveFluxEstimate;


namespace internal {


template< class DiffusionFactorImp, class DiffusiveFluxImp, class DiffusionTensorImp = void >
class DiffusiveFluxEstimateTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorImp >::value,
                "DiffusionFactorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxImp >::value,
                "DiffusiveFluxImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionTensorImp >::value,
                "DiffusionTensorImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_same< typename DiffusionFactorImp::EntityType,
                              typename DiffusiveFluxImp::EntityType >::value &&
                std::is_same< typename DiffusionFactorImp::EntityType,
                              typename DiffusionTensorImp::EntityType >::value,
                "EntityTypes have to agree!");
  static_assert(std::is_same< typename DiffusionFactorImp::DomainFieldType,
                              typename DiffusiveFluxImp::DomainFieldType >::value &&
                std::is_same< typename DiffusionFactorImp::DomainFieldType,
                              typename DiffusionTensorImp::DomainFieldType >::value,
                "DomainFieldTypes have to agree!");
  static_assert(DiffusionFactorImp::dimDomain == DiffusiveFluxImp::dimDomain &&
                DiffusionFactorImp::dimDomain == DiffusionTensorImp::dimDomain,
                "Dimensions of domains have to agree");
public:
  typedef DiffusiveFluxEstimate< DiffusionFactorImp, DiffusiveFluxImp, DiffusionTensorImp > derived_type;

  typedef DiffusionFactorImp                                          LocalizableDiffusionFactorType;
  typedef DiffusiveFluxImp                                            LocalizableDiffusiveFluxType;
  typedef DiffusionTensorImp                                          LocalizableDiffusionTensorType;

  typedef typename LocalizableDiffusionFactorType::LocalfunctionType  LocalDiffusionFactorType;
  typedef typename LocalizableDiffusiveFluxType::LocalfunctionType    LocalDiffusiveFluxType;
  typedef typename LocalizableDiffusionTensorType::LocalfunctionType  LocalDiffusionTensorType;
  typedef std::tuple< std::shared_ptr< LocalDiffusionFactorType >,
                      std::shared_ptr< LocalDiffusiveFluxType >,
                      std::shared_ptr< LocalDiffusionTensorType > >   LocalfunctionTupleType;
  typedef typename LocalizableDiffusionFactorType::EntityType         EntityType;
  typedef typename LocalizableDiffusionFactorType::DomainFieldType    DomainFieldType;
  static const unsigned int dimDomain = LocalizableDiffusionFactorType::dimDomain;
}; // class DiffusiveFluxEstimateTraits


template< class DiffusionImp, class DiffusiveFluxImp >
class DiffusiveFluxEstimateTraits< DiffusionImp, DiffusiveFluxImp, void >
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionImp >::value,
                "DiffusionImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxImp >::value,
                "DiffusiveFluxImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_same< typename DiffusionImp::EntityType, typename DiffusiveFluxImp::EntityType >::value,
                "EntityImps have to agree!");
  static_assert(std::is_same< typename DiffusionImp::DomainFieldType,
                              typename DiffusiveFluxImp::DomainFieldType >::value,
                "DomainFieldTypes have to agree!");
  static_assert(DiffusionImp::dimDomain == DiffusiveFluxImp::dimDomain, "Dimensions of domains have to agree");
public:
  typedef DiffusiveFluxEstimate< DiffusionImp, DiffusiveFluxImp >     derived_type;
  typedef DiffusionImp                                                LocalizableDiffusionType;
  typedef DiffusiveFluxImp                                            LocalizableDiffusiveFluxType;
  typedef typename LocalizableDiffusionType::LocalfunctionType        LocalDiffusionType;
  typedef typename LocalizableDiffusiveFluxType::LocalfunctionType    LocalDiffusiveFluxType;
  typedef std::tuple< std::shared_ptr< LocalDiffusionType >,
                      std::shared_ptr< LocalDiffusiveFluxType > >     LocalfunctionTupleType;
  typedef typename LocalizableDiffusionType::EntityType               EntityType;
  typedef typename LocalizableDiffusionType::DomainFieldType          DomainFieldType;
  static const unsigned int dimDomain = LocalizableDiffusionType::dimDomain;
}; // class DiffusiveFluxEstimateTraits


} // namespace internal


template< class DiffusionImp, class DiffusiveFluxImp >
class DiffusiveFluxEstimate< DiffusionImp, DiffusiveFluxImp, void >
  : public LocalEvaluation::Codim0Interface
         < internal::DiffusiveFluxEstimateTraits< DiffusionImp, DiffusiveFluxImp >, 2 >
{
public:
  typedef internal::DiffusiveFluxEstimateTraits< DiffusionImp, DiffusiveFluxImp > Traits;
  typedef typename Traits::LocalizableDiffusionType                               LocalizableDiffusionType;
  typedef typename Traits::LocalizableDiffusiveFluxType                           LocalizableDiffusiveFluxType;
  typedef typename Traits::LocalfunctionTupleType                                 LocalfunctionTupleType;
  typedef typename Traits::EntityType                                             EntityType;
  typedef typename Traits::DomainFieldType                                        DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  DiffusiveFluxEstimate(const LocalizableDiffusionType& diffusion, const LocalizableDiffusiveFluxType& diffusive_flux)
    : diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity), diffusive_flux_.local_function(entity));
  }

  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalfunctionTupleType localFuncs,
                      const Stuff::LocalfunctionSetInterface
                          < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
                      const Stuff::LocalfunctionSetInterface
                          < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase) const
  {
    const auto local_diffusion = std::get< 0 >(localFuncs);
    const auto local_diffusive_flux = std::get< 1 >(localFuncs);
    return redirect_order(*local_diffusion, *local_diffusive_flux, testBase, ansatzBase);
  } // ... order(...)

  template< class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalfunctionTupleType localFuncs,
                       const Stuff::LocalfunctionSetInterface
                           < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                       const Stuff::LocalfunctionSetInterface
                           < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base,
                       const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                       Dune::DynamicMatrix< R >& ret) const
  {
    const auto local_diffusion = std::get< 0 >(localFuncs);
    const auto local_diffusive_flux = std::get< 1 >(localFuncs);
    redirect_evaluate(*local_diffusion, *local_diffusive_flux, test_base, ansatz_base, local_point, ret);
  } // ... evaluate(...)

private:
  template< class R, int rLD, int rCLD, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  size_t redirect_order(const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& local_diffusion,
                        const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLDF, rCLDF >& local_diffusive_flux,
                        const Stuff::LocalfunctionSetInterface
                            < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                        const Stuff::LocalfunctionSetInterface
                            < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base) const
  {
    // TODO: there is no way to guess the order of local_diffusion^-1, so we take local_diffusion.order()
    return local_diffusion.order()
        + (std::max(local_diffusion.order() + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()))
        + (std::max(local_diffusion.order() + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()));
  } // ... redirect_order(...)

  template< class R, int rLD, int rCLD, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& /*local_diffusion*/,
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
  } // ... redirect_evaluate(...)

  template< class R >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& local_diffusion,
                         const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, dimDomain >& local_diffusive_flux,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& test_base,
                         const Stuff::LocalfunctionSetInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& ansatz_base,
                         const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                         Dune::DynamicMatrix< R >& ret) const
  {
    typedef typename Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, 1 >::JacobianRangeType JacobianRangeType;
    JacobianRangeType left_sum(0);
    JacobianRangeType right_sum(0);
    // evaluate local functions
    const auto diffusion_value = local_diffusion.evaluate(local_point);
    const R one_over_diffusion_value = 1.0 / diffusion_value;
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
  } // ... redirect_evaluate(...)

  const LocalizableDiffusionType& diffusion_;
  const LocalizableDiffusiveFluxType& diffusive_flux_;
}; // class DiffusiveFluxEstimate< ..., void >


template< class DiffusionFactorImp, class DiffusiveFluxImp, class DiffusionTensorImp >
class DiffusiveFluxEstimate
  : public LocalEvaluation::Codim0Interface< internal::DiffusiveFluxEstimateTraits< DiffusionFactorImp,
                                                                                    DiffusiveFluxImp,
                                                                                    DiffusionTensorImp >,
                                             2 >
{
public:
  typedef internal::DiffusiveFluxEstimateTraits< DiffusionFactorImp, DiffusiveFluxImp, DiffusionTensorImp > Traits;
  typedef typename Traits::LocalizableDiffusionFactorType                         LocalizableDiffusionFactorType;
  typedef typename Traits::LocalizableDiffusiveFluxType                           LocalizableDiffusiveFluxType;
  typedef typename Traits::LocalizableDiffusionTensorType                         LocalizableDiffusionTensorType;
  typedef typename Traits::LocalfunctionTupleType                                 LocalfunctionTupleType;
  typedef typename Traits::EntityType                                             EntityType;
  typedef typename Traits::DomainFieldType                                        DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  DiffusiveFluxEstimate(const LocalizableDiffusionFactorType& diffusion_factor,
                        const LocalizableDiffusionTensorType& diffusion_tensor,
                        const LocalizableDiffusiveFluxType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalfunctionTupleType localFuncs,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase) const
  {
    const auto local_diffusion_factor = std::get< 0 >(localFuncs);
    const auto local_diffusion_tensor = std::get< 1 >(localFuncs);
    const auto local_diffusive_flux = std::get< 2 >(localFuncs);
    return redirect_order(*local_diffusion_factor, *local_diffusion_tensor, *local_diffusive_flux, testBase,
                          ansatzBase);
  } // ... order(...)

  template< class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalfunctionTupleType localFuncs,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& test_base,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatz_base,
                const Dune::FieldVector< DomainFieldType, dimDomain >& local_point,
                Dune::DynamicMatrix< R >& ret) const
  {
    const auto local_diffusion_factor = std::get< 0 >(localFuncs);
    const auto local_diffusion_tensor = std::get< 1 >(localFuncs);
    const auto local_diffusive_flux = std::get< 2 >(localFuncs);
    redirect_evaluate(*local_diffusion_factor,
                      *local_diffusion_tensor,
                      *local_diffusive_flux,
                      test_base, ansatz_base, local_point, ret);
  } // ... evaluate(...)

private:
  template< class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  size_t redirect_order(const Stuff::LocalfunctionInterface
                            < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& local_diffusion_factor,
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
    const size_t local_diffusion_order = local_diffusion_factor.order() + local_diffusion_tensor.order();
    return local_diffusion_order
        + (std::max(local_diffusion_order + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()))
        + (std::max(local_diffusion_order + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()));
  } // ... redirect_order(...)

  template< class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, rLD, rCLD >& /*local_diffusion_factor*/,
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
  } // ... redirect_evaluate(...)

  template< class R >
  void redirect_evaluate(const Stuff::LocalfunctionInterface
                             < EntityType, DomainFieldType, dimDomain, R, 1 >& local_diffusion_factor,
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
    const auto diffusion_factor_value = local_diffusion_factor.evaluate(local_point);
    typedef Stuff::Common::FieldMatrix< R, dimDomain, dimDomain > TensorType;
    const TensorType diffusion_tensor_value = local_diffusion_tensor.evaluate(local_point);
    const TensorType diffusion_value = diffusion_tensor_value * diffusion_factor_value;
    // TODO: there is no documented way to assert that the inversion was successfull, so find one or check the matrix
    //       beforehand
    TensorType one_over_diffusion_value = diffusion_tensor_value;
    one_over_diffusion_value.invert();
    one_over_diffusion_value /= diffusion_factor_value;
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
        retRow[jj] = (one_over_diffusion_value * left_sum) * right_sum;
      }
    }
  } // ... redirect_evaluate(...)

  const LocalizableDiffusionFactorType& diffusion_factor_;
  const LocalizableDiffusionTensorType& diffusion_tensor_;
  const LocalizableDiffusiveFluxType& diffusive_flux_;
}; // class DiffusiveFluxEstimate


} // namespace ESV2007
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_ESV2007_HH
