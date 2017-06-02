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

#include <dune/xt/common/fmatrix.hh>
#include <dune/xt/grid/boundaryinfo/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/localizable-function.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/operators/integrals.hh>
#include <dune/gdt/local/integrands/lambda.hh>
#include <dune/gdt/operators/base.hh>
#include <dune/gdt/operators/oswaldinterpolation.hh>
#include <dune/gdt/spaces/dg/dune-fem-wrapper.hh>

namespace Dune {
namespace GDT {
namespace ESV2007 {

#if HAVE_DUNE_FEM


template <class GL, class InterpolationGridLayerType>
class NonconformityProduct
    : public LocalizableProductBase<GL,
                                    XT::Functions::LocalizableFunctionInterface<XT::Grid::extract_entity_t<GL>,
                                                                                typename GL::ctype,
                                                                                GL::dimension,
                                                                                double,
                                                                                1>>
{
  static_assert(XT::Grid::is_layer<GL>::value, "");
  static_assert(XT::Grid::is_layer<InterpolationGridLayerType>::value, "");
  typedef XT::Grid::extract_entity_t<GL> E;
  static_assert(std::is_same<XT::Grid::extract_entity_t<InterpolationGridLayerType>, E>::value, "");
  typedef typename GL::ctype D;
  static const constexpr size_t d = GL::dimension;
  typedef double R;

public:
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, 1> ScalarFunctionType;
  typedef XT::Functions::LocalizableFunctionInterface<E, D, d, R, d, d> TensorFunctionType;

private:
  typedef LocalizableProductBase<GL, ScalarFunctionType> BaseType;
  typedef NonconformityProduct<GL, InterpolationGridLayerType> ThisType;
  typedef LocalVolumeIntegralOperator<LocalLambdaBinaryVolumeIntegrand<E>,
                                      typename ScalarFunctionType::LocalfunctionType>
      LocalProductType;
  typedef DuneFemDgSpaceWrapper<InterpolationGridLayerType, 1, R, 1> DgSpaceType;
  typedef DiscreteFunction<DgSpaceType> DiscreteFunctionType;

public:
  using typename BaseType::GridLayerType;

  NonconformityProduct(GridLayerType grd_layr,
                       InterpolationGridLayerType interpolation_grid_layer,
                       const XT::Grid::BoundaryInfo<XT::Grid::extract_intersection_t<InterpolationGridLayerType>>&
                           interpolation_boundary_info,
                       const ScalarFunctionType& lambda,
                       const TensorFunctionType& kappa,
                       const ScalarFunctionType& u,
                       const ScalarFunctionType& v)
    : BaseType(grd_layr, u, v)
    , interpolation_grid_layer_(interpolation_grid_layer)
    , lambda_(lambda)
    , kappa_(kappa)
    , dg_space_(interpolation_grid_layer_)
    , interpolated_u_(dg_space_)
    , interpolated_v_(dg_space_)
    , local_product_(
          [&](const auto& local_u, const auto& local_v) {
            const auto& entity = local_u.entity();
            const auto local_lambda = lambda_.local_function(entity);
            const auto local_kappa = kappa_.local_function(entity);
            return local_lambda->order() + local_kappa->order()
                   + size_t(std::max(ssize_t(local_u.order()) - 1, ssize_t(0))
                            + std::max(ssize_t(local_v.order()) - 1, ssize_t(0)));
          },
          [&](const auto& local_u, const auto& local_v, const auto& local_point, auto& ret) {
            const auto& entity = local_u.entity();
            XT::Common::FieldMatrix<R, d, d> diffusion = kappa_.local_function(entity)->evaluate(local_point);
            diffusion *= lambda_.local_function(entity)->evaluate(local_point);
            const auto grad_u = local_u.jacobian(local_point).at(0)[0];
            const auto grad_interpolated_u = interpolated_u_.local_function(entity)->jacobian(local_point)[0];
            const auto grad_v = local_v.jacobian(local_point).at(0)[0];
            const auto grad_interpolated_v = interpolated_v_.local_function(entity)->jacobian(local_point)[0];
            ret[0][0] = (diffusion * (grad_u - grad_interpolated_u)) * (grad_v - grad_interpolated_v);
          })
  {
    OswaldInterpolationOperator<InterpolationGridLayerType> oswald_interpolation(interpolation_grid_layer_,
                                                                                 interpolation_boundary_info);
    oswald_interpolation.apply(this->range(), interpolated_u_);
    oswald_interpolation.apply(this->source(), interpolated_v_);
    this->append(local_product_);
  }

  NonconformityProduct(const ThisType&) = delete;
  NonconformityProduct(ThisType&&) = delete;

private:
  const InterpolationGridLayerType interpolation_grid_layer_;
  const ScalarFunctionType& lambda_;
  const TensorFunctionType& kappa_;
  const DgSpaceType dg_space_;
  DiscreteFunctionType interpolated_u_;
  DiscreteFunctionType interpolated_v_;
  const LocalProductType local_product_;
}; // class NonconformityProduct


#else // HAVE_DUNE_FEM


template <class GL, class InterpolationGridLayerType>
class NonconformityProduct
{
  static_assert(AlwaysFalse<GL>::value, "You are missing dune-fem!");
};


#endif // HAVE_DUNE_FEM

} // namespace ESV2007
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_OPERATORS_RS2017_HH
