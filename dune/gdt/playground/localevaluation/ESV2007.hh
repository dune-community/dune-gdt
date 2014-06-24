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


// forward, to be used int the traits
template< class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType = void >
class DiffusiveFluxEstimate;


namespace internal {


template< class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType = void >
class DiffusiveFluxEstimateTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionFactorType >::value,
                "DiffusionFactorType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxType >::value,
                "DiffusiveFluxType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionTensorType >::value,
                "DiffusionTensorType has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef DiffusiveFluxEstimate< DiffusionFactorType, DiffusiveFluxType, DiffusionTensorType > derived_type;
}; // class DiffusiveFluxEstimateTraits


template< class DiffusionType, class DiffusiveFluxType >
class DiffusiveFluxEstimateTraits< DiffusionType, DiffusiveFluxType, void >
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusionType >::value,
                "DiffusionType has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, DiffusiveFluxType >::value,
                "DiffusiveFluxType has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef DiffusiveFluxEstimate< DiffusionType, DiffusiveFluxType > derived_type;
}; // class DiffusiveFluxEstimateTraits


} // namespace internal


template< class DiffusionType, class DiffusiveFluxType >
class DiffusiveFluxEstimate< DiffusionType, DiffusiveFluxType, void >
  : public LocalEvaluation::Codim0Interface< internal::DiffusiveFluxEstimateTraits< DiffusionType, DiffusiveFluxType >,
                                             2 >
{
public:
  typedef internal::DiffusiveFluxEstimateTraits< DiffusionType, DiffusiveFluxType > Traits;

  DiffusiveFluxEstimate(const DiffusionType& diffusion, const DiffusiveFluxType& diffusive_flux)
    : diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename DiffusionType::LocalfunctionType     LocalDiffusionType;
    typedef typename DiffusiveFluxType::LocalfunctionType LocalDiffusiveFluxType;
  public:
    typedef std::tuple< std::shared_ptr< LocalDiffusionType >, std::shared_ptr< LocalDiffusiveFluxType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_.local_function(entity), diffusive_flux_.local_function(entity));
  }

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase)
  {
    const auto local_diffusion = std::get< 0 >(localFuncs);
    const auto local_diffusive_flux = std::get< 1 >(localFuncs);
    return order(*local_diffusion, *local_diffusive_flux, testBase, ansatzBase);
  } // ... order(...)

  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& local_diffusion,
                      const Stuff::LocalfunctionInterface< E, D, d, R, rLDF, rCLDF >& local_diffusive_flux,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base)
  {
    // TODO: there is no way to guess the order of local_diffusion^-1, so we take local_diffusion.order()
    return local_diffusion.order()
        + (std::max(local_diffusion.order() + std::max(ssize_t(test_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()))
        + (std::max(local_diffusion.order() + std::max(ssize_t(ansatz_base.order()) - 1, ssize_t(0)),
                    local_diffusive_flux.order()));
  } // ... order(...)

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base,
                       const Dune::FieldVector< D, d >& local_point,
                       Dune::DynamicMatrix< R >& ret)
  {
    const auto local_diffusion = std::get< 0 >(localFuncs);
    const auto local_diffusive_flux = std::get< 1 >(localFuncs);
    evaluate(*local_diffusion, *local_diffusive_flux, test_base, ansatz_base, local_point, ret);
  } // ... evaluate(...)

  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& /*local_diffusion*/,
                       const Stuff::LocalfunctionInterface< E, D, d, R, rLDF, rCLDF >& /*local_diffusive_flux*/,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& /*test_base*/,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& /*ansatz_base*/,
                       const Dune::FieldVector< D, d >& /*local_point*/,
                       Dune::DynamicMatrix< R >& /*ret*/)
  {
    static_assert(AlwaysFalse< R >::value, "Not implemented for these dimensions!");
  }

  template< class E, class D, int d, class R >
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1 >& local_diffusion,
                       const Stuff::LocalfunctionInterface< E, D, d, R, d >& local_diffusive_flux,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1 >& test_base,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, 1 >& ansatz_base,
                       const Dune::FieldVector< D, d >& local_point,
                       Dune::DynamicMatrix< R >& ret)
  {
    typedef typename Stuff::LocalfunctionSetInterface< E, D, d, R, 1 >::JacobianRangeType JacobianRangeType;
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
  } // ... evaluate(...)

private:
  const DiffusionType& diffusion_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class DiffusiveFluxEstimate< ..., void >


template< class DiffusionFactorType, class DiffusiveFluxType, class DiffusionTensorType >
class DiffusiveFluxEstimate
  : public LocalEvaluation::Codim0Interface< internal::DiffusiveFluxEstimateTraits< DiffusionFactorType,
                                                                                    DiffusiveFluxType,
                                                                                    DiffusionTensorType >,
                                             2 >
{
public:
  typedef internal::DiffusiveFluxEstimateTraits< DiffusionFactorType, DiffusiveFluxType, DiffusionTensorType > Traits;

  DiffusiveFluxEstimate(const DiffusionFactorType& diffusion_factor,
                        const DiffusionTensorType& diffusion_tensor,
                        const DiffusiveFluxType& diffusive_flux)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
  {}

  template< class EntityType >
  class LocalfunctionTuple
  {
    typedef typename DiffusionFactorType::LocalfunctionType LocalDiffusionFactorType;
    typedef typename DiffusionTensorType::LocalfunctionType LocalDiffusionTensorType;
    typedef typename DiffusiveFluxType::LocalfunctionType   LocalDiffusiveFluxType;
  public:
    typedef std::tuple< std::shared_ptr< LocalDiffusionFactorType >,
                        std::shared_ptr< LocalDiffusionTensorType >,
                        std::shared_ptr< LocalDiffusiveFluxType > > Type;
  };

  template< class EntityType >
  typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(diffusion_factor_.local_function(entity),
                           diffusion_tensor_.local_function(entity),
                           diffusive_flux_.local_function(entity));
  }

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static size_t order(const typename LocalfunctionTuple< E >::Type& localFuncs,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase)
  {
    const auto local_diffusion_factor = std::get< 0 >(localFuncs);
    const auto local_diffusion_tensor = std::get< 1 >(localFuncs);
    const auto local_diffusive_flux = std::get< 2 >(localFuncs);
    return order(*local_diffusion_factor, *local_diffusion_tensor, *local_diffusive_flux, testBase, ansatzBase);
  } // ... order(...)

  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static size_t order(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& local_diffusion_factor,
                      const Stuff::LocalfunctionInterface< E, D, d, R, rLDT, rCLDT >& local_diffusion_tensor,
                      const Stuff::LocalfunctionInterface< E, D, d, R, rLDF, rCLDF >& local_diffusive_flux,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                      const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base)
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

  template< class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  static void evaluate(const typename LocalfunctionTuple< E >::Type& localFuncs,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& test_base,
                       const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatz_base,
                       const Dune::FieldVector< D, d >& local_point,
                       Dune::DynamicMatrix< R >& ret)
  {
    const auto local_diffusion_factor = std::get< 0 >(localFuncs);
    const auto local_diffusion_tensor = std::get< 1 >(localFuncs);
    const auto local_diffusive_flux = std::get< 2 >(localFuncs);
    evaluate(*local_diffusion_factor,
             *local_diffusion_tensor,
             *local_diffusive_flux,
             test_base, ansatz_base, local_point, ret);
  } // ... evaluate(...)

  template< class E, class D, int d, class R, int rLD, int rCLD, int rLDT, int rCLDT, int rLDF, int rCLDF, int rT, int rCT, int rA, int rCA >
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, rLD, rCLD >& /*local_diffusion_factor*/,
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
  static void evaluate(const Stuff::LocalfunctionInterface< E, D, d, R, 1 >& local_diffusion_factor,
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
    const auto diffusion_factor_value = local_diffusion_factor.evaluate(local_point);
    typedef Stuff::Common::FieldMatrix< R, d, d > TensorType;
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
  } // ... evaluate(...)

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
}; // class DiffusiveFluxEstimate


} // namespace ESV2007
} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALEVALUATION_ESV2007_HH
