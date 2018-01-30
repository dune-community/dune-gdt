// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PLAYGROUND_OPERATORS_OS2015_HH
#define DUNE_GDT_PLAYGROUND_OPERATORS_OS2015_HH

#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/la/container/eigen.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/entity.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/derived.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/local/integrands/lambda.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/fluxreconstruction.hh>
#include <dune/gdt/spaces/rt/default.hh>

namespace Dune {
namespace GDT {
namespace OS2015 {
namespace internal {


template <class ProductGridLayer, class ReconstructionGridLayer>
class ResidualProductBase
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
  typedef ResidualProductBase<ProductGridLayer, ReconstructionGridLayer> ThisType;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef RaviartThomasSpace<ReconstructionGridLayer, 0, R> RtSpaceType;
  typedef DiscreteFunction<RtSpaceType> FluxReconstructionType;
  typedef XT::Functions::DivergenceFunction<FluxReconstructionType> DivergenceOfFluxReconstructionType;
  typedef typename ScalarFunctionType::DifferenceType DifferenceType;

public:
  ResidualProductBase(ReconstructionGridLayer reconstruction_grid_layer,
                      const ScalarFunctionType& lambda,
                      const TensorFunctionType& kappa,
                      const ScalarFunctionType& f,
                      const ScalarFunctionType& u,
                      const ScalarFunctionType& v)
    : f_(f)
    , rt_space_(reconstruction_grid_layer)
    , reconstructed_u_(rt_space_)
    , reconstructed_v_(rt_space_)
    , divergence_of_reconstructed_u_(reconstructed_u_)
    , divergence_of_reconstructed_v_(reconstructed_v_)
    , f_minus_divergence_of_reconstructed_u_(f_ - divergence_of_reconstructed_u_)
    , f_minus_divergence_of_reconstructed_v_(f_ - divergence_of_reconstructed_v_)
  {
    DiffusiveFluxReconstructionOperator<ReconstructionGridLayer,
                                        ScalarFunctionType,
                                        TensorFunctionType,
                                        LocalEllipticIpdgIntegrands::Method::swipdg_affine_factor>
        flux_reconstruction(reconstruction_grid_layer, lambda, kappa);
    flux_reconstruction.apply(u, reconstructed_u_);
    flux_reconstruction.apply(v, reconstructed_v_);
  }

  ResidualProductBase(const ThisType&) = delete;
  ResidualProductBase(ThisType&&) = delete;

protected:
  const ScalarFunctionType& f_;
  const RtSpaceType rt_space_;
  FluxReconstructionType reconstructed_u_;
  FluxReconstructionType reconstructed_v_;
  const DivergenceOfFluxReconstructionType divergence_of_reconstructed_u_;
  const DivergenceOfFluxReconstructionType divergence_of_reconstructed_v_;
  const DifferenceType f_minus_divergence_of_reconstructed_u_;
  const DifferenceType f_minus_divergence_of_reconstructed_v_;
}; // class ResidualProductBase


} // namespace internal


template <class ProductGridLayer, class ReconstructionGridLayer>
class ResidualProduct
    : internal::ResidualProductBase<ProductGridLayer, ReconstructionGridLayer>,
      public LocalizableProductBase<ProductGridLayer,
                                    XT::Functions::
                                        LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                     typename ProductGridLayer::ctype,
                                                                     ProductGridLayer::dimension,
                                                                     double,
                                                                     1>>
{
  typedef internal::ResidualProductBase<ProductGridLayer, ReconstructionGridLayer> ResidualProductBaseType;
  typedef LocalizableProductBase<ProductGridLayer,
                                 XT::Functions::
                                     LocalizableFunctionInterface<XT::Grid::extract_entity_t<ProductGridLayer>,
                                                                  typename ProductGridLayer::ctype,
                                                                  ProductGridLayer::dimension,
                                                                  double,
                                                                  1>>
      LocalizableProductBaseType;

public:
  using typename ResidualProductBaseType::ScalarFunctionType;
  using typename ResidualProductBaseType::TensorFunctionType;

private:
  using typename ResidualProductBaseType::E;
  using typename ResidualProductBaseType::D;
  using ResidualProductBaseType::d;
  using typename ResidualProductBaseType::R;
  typedef LocalVolumeIntegralOperator<LocalLambdaBinaryVolumeIntegrand<E>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;

public:
  ResidualProduct(ProductGridLayer product_grid_layer,
                  ReconstructionGridLayer reconstruction_grid_layer,
                  const ScalarFunctionType& lambda,
                  const ScalarFunctionType& lambda_hat,
                  const TensorFunctionType& kappa,
                  const ScalarFunctionType& f,
                  const ScalarFunctionType& u,
                  const ScalarFunctionType& v,
                  const double& poincare_constant = 1.0 / (M_PIl * M_PIl),
                  const size_t over_integrate = 2)
    : ResidualProductBaseType(reconstruction_grid_layer, lambda, kappa, f, u, v)
    , LocalizableProductBaseType(product_grid_layer,
                                 this->f_minus_divergence_of_reconstructed_u_,
                                 this->f_minus_divergence_of_reconstructed_v_)
    , lambda_(lambda_hat)
    , kappa_(kappa)
    , poincare_constant_(poincare_constant)
    , over_integrate_(over_integrate)
    , local_product_(
          // the order lambda
          [&](const auto& local_f_minus_divergence_of_reconstructed_u,
              const auto& local_f_minus_divergence_of_reconstructed_v) {
            return local_f_minus_divergence_of_reconstructed_u.order()
                   + local_f_minus_divergence_of_reconstructed_v.order() + over_integrate_;
          },
          // the evaluate lambda
          [&](const auto& local_f_minus_divergence_of_reconstructed_u,
              const auto& local_f_minus_divergence_of_reconstructed_v,
              const auto& local_point,
              auto& ret) {
            ret[0][0] = local_f_minus_divergence_of_reconstructed_u.evaluate(local_point).at(0)[0]
                        * local_f_minus_divergence_of_reconstructed_v.evaluate(local_point).at(0)[0];
          })
    , min_diffusion_ev_(std::numeric_limits<R>::max())
    , subdomain_vertices_()
    , finalized_(false)
  {
    this->append(local_product_);

    // the functor to collect all grid vertices
    this->append([&](const E& entity) {
      for (size_t cc = 0; cc < entity.subEntities(d); ++cc)
        subdomain_vertices_.emplace_back(entity.template subEntity<d>(cc).geometry().center());
    });

    // the functor to compute the min eigenvalue of the diffusion
    this->append([&](const E& entity) {
      const auto local_lambda = lambda_.local_function(entity);
      const auto local_kappa = kappa_.local_function(entity);
      // To find the minimum of a function we evaluate it
      // * in all quadrature points of a quadrature which would integrate such a function exactly
      for (const auto& quadrature_point :
           QuadratureRules<D, d>::rule(entity.type(), local_lambda->order() + local_kappa->order() + over_integrate_)) {
        const auto xx = quadrature_point.position();
        auto diffusion = local_kappa->evaluate(xx);
        diffusion *= local_lambda->evaluate(xx);
        min_diffusion_ev_ = std::min(
            min_diffusion_ev_,
            XT::LA::make_eigen_solver(diffusion, {"assert_positive_eigenvalues", 1e-15}).min_eigenvalues(1).at(0));
      }
      // * and in the corners of the gigen entity.
      const auto& reference_element = ReferenceElements<D, d>::general(entity.type());
      for (int ii = 0; ii < reference_element.size(d); ++ii) {
        const auto xx = reference_element.position(ii, d);
        auto diffusion = local_kappa->evaluate(xx);
        diffusion *= local_lambda->evaluate(xx);
        min_diffusion_ev_ = std::min(
            min_diffusion_ev_,
            XT::LA::make_eigen_solver(diffusion, {"assert_positive_eigenvalues", 1e-15}).min_eigenvalues(1).at(0));
      }
    });
  } // ResidualProduct(...)

  R apply2()
  {
    if (!finalized_) {
      this->walk();
      R subdomain_h = std::numeric_limits<R>::min();
      for (size_t ii = 0; ii < subdomain_vertices_.size(); ++ii)
        for (size_t jj = ii + 1; jj < subdomain_vertices_.size(); ++jj)
          subdomain_h = std::max(subdomain_h, (subdomain_vertices_[ii] - subdomain_vertices_[jj]).two_norm());
      this->result_ *= (poincare_constant_ / min_diffusion_ev_) * subdomain_h * subdomain_h;
    }
    return this->result_;
  }

private:
  const ScalarFunctionType& lambda_;
  const TensorFunctionType& kappa_;
  const double poincare_constant_;
  const size_t over_integrate_;
  const LocalProductType local_product_;
  R min_diffusion_ev_;
  std::vector<FieldVector<D, d>> subdomain_vertices_;
  bool finalized_;
}; // class ResidualProduct


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
  typedef RaviartThomasSpace<ReconstructionGridLayer, 0, R> RtSpaceType;
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


} // namespace OS2015
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_OPERATORS_OS2015_HH
