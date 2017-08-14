// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PLAYGROUND_OPERATORS_RS2017_HH
#define DUNE_GDT_PLAYGROUND_OPERATORS_RS2017_HH

#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/derived.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/local/integrands/lambda.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/fluxreconstruction.hh>
#include <dune/gdt/spaces/rt/dune-pdelab-wrapper.hh>

namespace Dune {
namespace GDT {
namespace RS2017 {
namespace DiffusiveFluxIndicator {
#if HAVE_DUNE_PDELAB

#if 0
namespace ResidualIndicator {


template <class ProductGridLayer>
class FtimesFproduct
    : public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  static_assert(XT::Grid::is_layer<ProductGridLayer>::value, "");
  typedef XT::Grid::extract_entity_t<ProductGridLayer> E;
  typedef typename ProductGridLayer::ctype D;
  static const constexpr size_t d = ProductGridLayer::dimension;
  typedef double R;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;

private:
  typedef LocalizableProductBase<ProductGridLayer, ScalarFunctionType> BaseType;
  typedef FtimesFproduct<ProductGridLayer> ThisType;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<ScalarFunctionType>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  using typename BaseType::GridLayerType;

  FtimesFproduct(GridLayerType grd_layr, const ScalarFunctionType& f, const size_t over_integrate = 2)
    : BaseType(grd_layr, f, f)
    , one_(1.)
    , local_product_(over_integrate, one_)
  {
    this->append(local_product_);
  }

  FtimesFproduct(const ThisType&) = delete;
  FtimesFproduct(ThisType&&) = delete;

private:
  const XT::Functions::ConstantFunction<E, D, d, R, 1> one_;
  const LocalProductType local_product_;
}; // class FtimesFproduct


namespace internal {


template <class ProductGridLayer, class ReconstructionGridLayer>
class FtimesVproductBase
{
  static_assert(XT::Grid::is_layer<ProductGridLayer>::value, "");
  static_assert(XT::Grid::is_layer<ReconstructionGridLayer>::value, "");

protected:
  typedef XT::Grid::extract_entity_t<ProductGridLayer> E;
  typedef typename ProductGridLayer::ctype D;
  static const constexpr size_t d = ProductGridLayer::dimension;
  typedef double R;

private:
  static_assert(std::is_same<XT::Grid::extract_entity_t<ReconstructionGridLayer>, E>::value, "");
  typedef FtimesVproductBase<ProductGridLayer, ReconstructionGridLayer> ThisType;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef DunePdelabRtSpaceWrapper<ReconstructionGridLayer, 0, R, d> RtSpaceType;
  typedef DiscreteFunction<RtSpaceType> FluxReconstructionType;
  typedef XT::Functions::DivergenceFunction<FluxReconstructionType> DivergenceOfFluxReconstructionType;
  typedef typename ScalarFunctionType::DifferenceType DifferenceType;

public:
  FtimesVproductBase(ReconstructionGridLayer reconstruction_grid_layer,
                     const ScalarFunctionType& lambda,
                     const TensorFunctionType& kappa,
                     const ScalarFunctionType& f,
                     const ScalarFunctionType& v)
    : f_(f)
    , rt_space_(reconstruction_grid_layer)
    , reconstructed_v_(rt_space_)
    , divergence_of_reconstructed_v_(reconstructed_v_)
  {
    DiffusiveFluxReconstructionOperator<ReconstructionGridLayer,
                                        ScalarFunctionType,
                                        TensorFunctionType,
                                        LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        flux_reconstruction(reconstruction_grid_layer, lambda, kappa);
    flux_reconstruction.apply(v, reconstructed_v_);
  }

  FtimesVproductBase(const ThisType&) = delete;
  FtimesVproductBase(ThisType&&) = delete;

protected:
  const ScalarFunctionType& f_;
  const RtSpaceType rt_space_;
  FluxReconstructionType reconstructed_v_;
  const DivergenceOfFluxReconstructionType divergence_of_reconstructed_v_;
}; // class FtimesVproductBase


template <class ProductGridLayer, class ReconstructionGridLayer>
class UtimesVproductBase
{
  static_assert(XT::Grid::is_layer<ProductGridLayer>::value, "");
  static_assert(XT::Grid::is_layer<ReconstructionGridLayer>::value, "");

protected:
  typedef XT::Grid::extract_entity_t<ProductGridLayer> E;
  typedef typename ProductGridLayer::ctype D;
  static const constexpr size_t d = ProductGridLayer::dimension;
  typedef double R;

private:
  static_assert(std::is_same<XT::Grid::extract_entity_t<ReconstructionGridLayer>, E>::value, "");
  typedef UtimesVproductBase<ProductGridLayer, ReconstructionGridLayer> ThisType;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef DunePdelabRtSpaceWrapper<ReconstructionGridLayer, 0, R, d> RtSpaceType;
  typedef DiscreteFunction<RtSpaceType> FluxReconstructionType;
  typedef XT::Functions::DivergenceFunction<FluxReconstructionType> DivergenceOfFluxReconstructionType;
  typedef typename ScalarFunctionType::DifferenceType DifferenceType;

public:
  UtimesVproductBase(ReconstructionGridLayer reconstruction_grid_layer,
                     const ScalarFunctionType& lambda,
                     const TensorFunctionType& kappa,
                     const ScalarFunctionType& u,
                     const ScalarFunctionType& v)
    : rt_space_(reconstruction_grid_layer)
    , reconstructed_u_(rt_space_)
    , reconstructed_v_(rt_space_)
    , divergence_of_reconstructed_u_(reconstructed_u_)
    , divergence_of_reconstructed_v_(reconstructed_v_)
  {
    DiffusiveFluxReconstructionOperator<ReconstructionGridLayer,
                                        ScalarFunctionType,
                                        TensorFunctionType,
                                        LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        flux_reconstruction(reconstruction_grid_layer, lambda, kappa);
    flux_reconstruction.apply(u, reconstructed_u_);
    flux_reconstruction.apply(v, reconstructed_v_);
  }

  UtimesVproductBase(const ThisType&) = delete;
  UtimesVproductBase(ThisType&&) = delete;

protected:
  const RtSpaceType rt_space_;
  FluxReconstructionType reconstructed_u_;
  FluxReconstructionType reconstructed_v_;
  const DivergenceOfFluxReconstructionType divergence_of_reconstructed_u_;
  const DivergenceOfFluxReconstructionType divergence_of_reconstructed_v_;
}; // class UtimesVproductBase


} // namespace internal


template <class ProductGridLayer, class ReconstructionGridLayer>
class FtimesVproduct
    : internal::FtimesVproductBase<ProductGridLayer, ReconstructionGridLayer>,
      public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  typedef internal::FtimesVproductBase<ProductGridLayer, ReconstructionGridLayer> FtimesVproductBaseType;
  typedef LocalizableProductBase<ProductGridLayer,
                                 XT::Functions::
                                     LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                  typename ProductGridLayer::ctype,
                                                                  ProductGridLayer::dimension,
                                                                  double,
                                                                  1>>
      LocalizableProductBaseType;

public:
  using typename FtimesVproductBaseType::ScalarFunctionType;
  using typename FtimesVproductBaseType::TensorFunctionType;

private:
  using typename FtimesVproductBaseType::E;
  using typename FtimesVproductBaseType::D;
  using FtimesVproductBaseType::d;
  using typename FtimesVproductBaseType::R;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<ScalarFunctionType>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  FtimesVproduct(ProductGridLayer product_grid_layer,
                 ReconstructionGridLayer reconstruction_grid_layer,
                 const ScalarFunctionType& lambda,
                 const TensorFunctionType& kappa,
                 const ScalarFunctionType& f,
                 const ScalarFunctionType& v,
                 const size_t over_integrate = 2)
    : FtimesVproductBaseType(reconstruction_grid_layer, lambda, kappa, f, v)
    , LocalizableProductBaseType(product_grid_layer, this->f_, this->divergence_of_reconstructed_v_)
    , one_(1.)
    , local_product_(over_integrate, one_)
  {
    this->append(local_product_);
  }

private:
  const XT::Functions::ConstantFunction<E, D, d, R, 1> one_;
  const LocalProductType local_product_;
}; // class FtimesVproduct


template <class ProductGridLayer, class ReconstructionGridLayer>
class UtimesVproduct
    : internal::UtimesVproductBase<ProductGridLayer, ReconstructionGridLayer>,
      public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  typedef internal::UtimesVproductBase<ProductGridLayer, ReconstructionGridLayer> UtimesVproductBaseType;
  typedef LocalizableProductBase<ProductGridLayer,
                                 XT::Functions::
                                     LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                  typename ProductGridLayer::ctype,
                                                                  ProductGridLayer::dimension,
                                                                  double,
                                                                  1>>
      LocalizableProductBaseType;

public:
  using typename UtimesVproductBaseType::ScalarFunctionType;
  using typename UtimesVproductBaseType::TensorFunctionType;

private:
  using typename UtimesVproductBaseType::E;
  using typename UtimesVproductBaseType::D;
  using UtimesVproductBaseType::d;
  using typename UtimesVproductBaseType::R;
  typedef LocalVolumeIntegralOperator<LocalProductIntegrand<ScalarFunctionType>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  UtimesVproduct(ProductGridLayer product_grid_layer,
                 ReconstructionGridLayer reconstruction_grid_layer,
                 const ScalarFunctionType& lambda,
                 const TensorFunctionType& kappa,
                 const ScalarFunctionType& u,
                 const ScalarFunctionType& v,
                 const size_t over_integrate = 2)
    : UtimesVproductBaseType(reconstruction_grid_layer, lambda, kappa, u, v)
    , LocalizableProductBaseType(
          product_grid_layer, this->divergence_of_reconstructed_u_, this->divergence_of_reconstructed_v_)
    , one_(1.)
    , local_product_(over_integrate, one_)
  {
    this->append(local_product_);
  }

private:
  const XT::Functions::ConstantFunction<E, D, d, R, 1> one_;
  const LocalProductType local_product_;
}; // class UtimesVproduct


#else // HAVE_EIGEN


template <class ProductGridLayer, class ReconstructionGridLayer>
class FtimesVproduct
{
  static_assert(AlwaysFalse<ProductGridLayer>::value, "You are missing eigen!");
};


template <class ProductGridLayer, class ReconstructionGridLayer>
class UtimesVproduct
{
  static_assert(AlwaysFalse<ProductGridLayer>::value, "You are missing eigen!");
};


#endif // HAVE_EIGEN

} // namespace ResidualIndicator
#endif // 0


template <class ProductGridLayer, class ReconstructionGridLayer>
class AbProduct
    : public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  static_assert(XT::Grid::is_layer<ProductGridLayer>::value, "");
  static_assert(XT::Grid::is_layer<ReconstructionGridLayer>::value, "");
  typedef XT::Grid::extract_entity_t<ProductGridLayer> E;
  static_assert(std::is_same<XT::Grid::extract_entity_t<ReconstructionGridLayer>, E>::value, "");
  typedef typename ProductGridLayer::ctype D;
  static const constexpr size_t d = ProductGridLayer::dimension;
  typedef double R;
  typedef AbProduct<ProductGridLayer, ReconstructionGridLayer> ThisType;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef DunePdelabRtSpaceWrapper<ReconstructionGridLayer, 0, R, d> RtSpaceType;
  typedef DiscreteFunction<RtSpaceType> FluxReconstructionType;
  typedef LocalizableProductBase<ProductGridLayer, XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>> BaseType;
  typedef LocalVolumeIntegralOperator<LocalLambdaBinaryVolumeIntegrand<E>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  AbProduct(ProductGridLayer product_grid_layer,
            ReconstructionGridLayer reconstruction_grid_layer,
            const ScalarFunctionType& lambda,
            const ScalarFunctionType& lambda_hat,
            const TensorFunctionType& kappa,
            const ScalarFunctionType& u,
            const ScalarFunctionType& v,
            const size_t over_integrate = 2)
    : BaseType(product_grid_layer, u, v)
    , lambda_(lambda)
    , lambda_hat_(lambda_hat)
    , kappa_(kappa)
    , rt_space_(reconstruction_grid_layer)
    , reconstructed_u_(rt_space_)
    , reconstructed_v_(rt_space_)
    , over_integrate_(over_integrate)
    , local_product_(
          // the order lambda
          [&](const auto& local_u, const auto& local_v) {
            const auto& entity = local_u.entity();
            const size_t diffusion_order =
                lambda_.local_function(entity)->order() + kappa_.local_function(entity)->order();
            return 3 * diffusion_order + size_t(std::max(ssize_t(local_u.order()) - 1, ssize_t(0)))
                   + size_t(std::max(ssize_t(local_v.order()) - 1, ssize_t(0))) + over_integrate_;
          },
          // the evaluate lambda
          [&](const auto& local_u, const auto& local_v, const auto& local_point, auto& ret) {
            const auto& entity = local_u.entity();
            XT::Common::FieldMatrix<R, d, d> diffusion = kappa_.local_function(entity)->evaluate(local_point);
            XT::Common::FieldMatrix<R, d, d> one_over_diffusion_hat = diffusion;
            diffusion *= lambda_.local_function(entity)->evaluate(local_point);
            one_over_diffusion_hat *= lambda_hat_.local_function(entity)->evaluate(local_point);
            one_over_diffusion_hat.invert(); // there is no documented way to assert that the inversion was successfull
            const auto grad_u = local_u.jacobian(local_point).at(0)[0];
            const auto grad_v = local_v.jacobian(local_point).at(0)[0];
            const auto val_reconstructed_u = reconstructed_u_.local_function(entity)->evaluate(local_point);
            const auto val_reconstructed_v = reconstructed_v_.local_function(entity)->evaluate(local_point);
            ret[0][0] = (one_over_diffusion_hat * ((diffusion * grad_u) + val_reconstructed_u)) // clang-format off
                                                * ((diffusion * grad_v) + val_reconstructed_v); // clang-format on
          })
  {
    DiffusiveFluxReconstructionOperator<ReconstructionGridLayer,
                                        ScalarFunctionType,
                                        TensorFunctionType,
                                        LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        flux_reconstruction(reconstruction_grid_layer, lambda, kappa);
    flux_reconstruction.apply(u, reconstructed_u_);
    flux_reconstruction.apply(v, reconstructed_v_);
    this->append(local_product_);
  }

  AbProduct(const ThisType&) = delete;
  AbProduct(ThisType&&) = delete;

private:
  const ScalarFunctionType& lambda_;
  const ScalarFunctionType& lambda_hat_;
  const TensorFunctionType& kappa_;
  const RtSpaceType rt_space_;
  FluxReconstructionType reconstructed_u_;
  FluxReconstructionType reconstructed_v_;
  const size_t over_integrate_;
  const LocalProductType local_product_;
}; // class AbProduct


template <class ProductGridLayer, class ReconstructionGridLayer>
class DiffusiveFluxProduct
    : public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  static_assert(XT::Grid::is_layer<ProductGridLayer>::value, "");
  static_assert(XT::Grid::is_layer<ReconstructionGridLayer>::value, "");
  typedef XT::Grid::extract_entity_t<ProductGridLayer> E;
  static_assert(std::is_same<XT::Grid::extract_entity_t<ReconstructionGridLayer>, E>::value, "");
  typedef typename ProductGridLayer::ctype D;
  static const constexpr size_t d = ProductGridLayer::dimension;
  typedef double R;
  typedef DiffusiveFluxProduct<ProductGridLayer, ReconstructionGridLayer> ThisType;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef DunePdelabRtSpaceWrapper<ReconstructionGridLayer, 0, R, d> RtSpaceType;
  typedef DiscreteFunction<RtSpaceType> FluxReconstructionType;
  typedef LocalizableProductBase<ProductGridLayer, XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1>> BaseType;
  typedef LocalVolumeIntegralOperator<LocalLambdaBinaryVolumeIntegrand<E>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  DiffusiveFluxProduct(ProductGridLayer product_grid_layer,
                       ReconstructionGridLayer reconstruction_grid_layer,
                       const ScalarFunctionType& lambda,
                       const ScalarFunctionType& lambda_hat,
                       const TensorFunctionType& kappa,
                       const ScalarFunctionType& u,
                       const ScalarFunctionType& v,
                       const size_t over_integrate = 2)
    : BaseType(product_grid_layer, u, v)
    , lambda_(lambda)
    , lambda_hat_(lambda_hat)
    , kappa_(kappa)
    , rt_space_(reconstruction_grid_layer)
    , reconstructed_u_(rt_space_)
    , reconstructed_v_(rt_space_)
    , over_integrate_(over_integrate)
    , local_product_(
          // the order lambda
          [&](const auto& local_u, const auto& local_v) {
            const auto& entity = local_u.entity();
            const size_t diffusion_order =
                lambda_.local_function(entity)->order() + kappa_.local_function(entity)->order();
            return 3 * diffusion_order + size_t(std::max(ssize_t(local_u.order()) - 1, ssize_t(0)))
                   + size_t(std::max(ssize_t(local_v.order()) - 1, ssize_t(0))) + over_integrate_;
          },
          // the evaluate lambda
          [&](const auto& local_u, const auto& local_v, const auto& local_point, auto& ret) {
            const auto& entity = local_u.entity();
            XT::Common::FieldMatrix<R, d, d> diffusion = kappa_.local_function(entity)->evaluate(local_point);
            XT::Common::FieldMatrix<R, d, d> one_over_diffusion_hat = diffusion;
            diffusion *= lambda_.local_function(entity)->evaluate(local_point);
            one_over_diffusion_hat *= lambda_hat_.local_function(entity)->evaluate(local_point);
            one_over_diffusion_hat.invert(); // there is no documented way to assert that the inversion was successfull
            const auto grad_u = local_u.jacobian(local_point).at(0)[0];
            const auto grad_v = local_v.jacobian(local_point).at(0)[0];
            const auto val_reconstructed_u = reconstructed_u_.local_function(entity)->evaluate(local_point);
            const auto val_reconstructed_v = reconstructed_v_.local_function(entity)->evaluate(local_point);
            ret[0][0] = (one_over_diffusion_hat * ((diffusion * grad_u) + val_reconstructed_u)) // clang-format off
                                                * ((diffusion * grad_v) + val_reconstructed_v); // clang-format on
          })
  {
    DiffusiveFluxReconstructionOperator<ReconstructionGridLayer,
                                        ScalarFunctionType,
                                        TensorFunctionType,
                                        LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        flux_reconstruction(reconstruction_grid_layer, lambda, kappa);
    flux_reconstruction.apply(u, reconstructed_u_);
    flux_reconstruction.apply(v, reconstructed_v_);
    this->append(local_product_);
  }

  DiffusiveFluxProduct(const ThisType&) = delete;
  DiffusiveFluxProduct(ThisType&&) = delete;

private:
  const ScalarFunctionType& lambda_;
  const ScalarFunctionType& lambda_hat_;
  const TensorFunctionType& kappa_;
  const RtSpaceType rt_space_;
  FluxReconstructionType reconstructed_u_;
  FluxReconstructionType reconstructed_v_;
  const size_t over_integrate_;
  const LocalProductType local_product_;
}; // class DiffusiveFluxProduct


#else // HAVE_DUNE_PDELAB


template <class ProductGridLayer, class ReconstructionGridLayer>
class DiffusiveFluxProduct
{
  static_assert(AlwaysFalse<ProductGridLayer>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB
} // namespace DiffusiveFluxIndicator
} // namespace RS2017
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_OPERATORS_RS2017_HH
