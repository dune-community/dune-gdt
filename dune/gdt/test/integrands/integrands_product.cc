// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/local/integrands/product.hh>

#include <dune/gdt/test/integrands/integrands.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct ProductIntegrandTest : public IntegrandTest<G>
{
  using BaseType = IntegrandTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::VectorJacobianType;
  using ScalarIntegrandType = LocalElementProductIntegrand<E, 1>;
  using VectorIntegrandType = LocalElementProductIntegrand<E, d>;

  virtual void is_constructable() override final
  {
    ScalarIntegrandType scalar_integrand1;
    ScalarIntegrandType scalar_integrand2(1.);
    const XT::Functions::GenericGridFunction<E, 1> scalar_grid_function(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    ScalarIntegrandType scalar_integrand3(scalar_grid_function);
    const XT::Functions::GenericFunction<d, 1> scalar_function(
        2, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    ScalarIntegrandType scalar_integrand4(scalar_function);
    DUNE_UNUSED_PARAMETER(scalar_integrand1);
    DUNE_UNUSED_PARAMETER(scalar_integrand2);
    DUNE_UNUSED_PARAMETER(scalar_integrand3);
    DUNE_UNUSED_PARAMETER(scalar_integrand4);
    VectorIntegrandType vector_integrand1;
    VectorIntegrandType vector_integrand2(1.);
    VectorIntegrandType vector_integrand3(scalar_grid_function);
    VectorIntegrandType vector_integrand4(scalar_function);
    const XT::Functions::GenericGridFunction<E, 2, 2> matrix_grid_function(
        1,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return VectorJacobianType{{x[0], x[1]}, {1., 2.}};
        });
    VectorIntegrandType vector_integrand5(matrix_grid_function);
    DUNE_UNUSED_PARAMETER(vector_integrand1);
    DUNE_UNUSED_PARAMETER(vector_integrand2);
    DUNE_UNUSED_PARAMETER(vector_integrand3);
    DUNE_UNUSED_PARAMETER(vector_integrand4);
    DUNE_UNUSED_PARAMETER(vector_integrand5);
  }

  virtual void evaluates_correctly_for_scalar_bases()
  {
    const XT::Functions::GenericGridFunction<E, 1> scalar_inducing_function(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    ScalarIntegrandType scalar_integrand(scalar_inducing_function);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    scalar_integrand.bind(element);
    const auto integrand_order = scalar_integrand.order(*scalar_test_, *scalar_ansatz_);
    EXPECT_EQ(9, integrand_order);
    DynamicMatrix<D> result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      scalar_integrand.evaluate(*scalar_test_, *scalar_ansatz_, x, result);
      DynamicMatrix<D> expected_result{{std::pow(x[0] * x[1], 2), std::pow(x[0] * x[1], 3)},
                                       {std::pow(x[0], 3) * std::pow(x[1], 4), std::pow(x[0], 4) * std::pow(x[1], 5)}};
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < 2; ++jj)
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
    }
  }

  virtual void evaluates_correctly_for_vector_bases()
  {
    const XT::Functions::GenericGridFunction<E, 2, 2> inducing_function(
        1,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return VectorJacobianType{{x[0], x[1]}, {1., 2.}};
        });
    VectorIntegrandType integrand(inducing_function);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    integrand.bind(element);
    const auto integrand_order = integrand.order(*vector_test_, *vector_ansatz_);
    EXPECT_EQ(5, integrand_order);
    DynamicMatrix<D> result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      integrand.evaluate(*vector_test_, *vector_ansatz_, x, result);
      DynamicMatrix<D> expected_result{{std::pow(x[0], 2) + 2 * x[0] * x[1] + 5 * x[1],
                                        std::pow(x[0], 3) + 2 * std::pow(x[0], 2) * x[1] + 5 * std::pow(x[1], 2)},
                                       {std::pow(x[0], 3) * x[1] + std::pow(x[0], 2) * x[1]
                                            + 2 * x[0] * std::pow(x[1], 2) + 2 * x[0] * x[1] + 2 * std::pow(x[1], 2),
                                        std::pow(x[0], 4) * x[1] + std::pow(x[0], 3) * x[1]
                                            + std::pow(x[0], 2) * std::pow(x[1], 2) + x[0] * std::pow(x[1], 3)
                                            + +2 * x[0] * std::pow(x[1], 2) + 2 * std::pow(x[1], 3)}};
      for (size_t ii = 0; ii < 2; ++ii)
        for (size_t jj = 0; jj < 2; ++jj)
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
    }
  }

  using BaseType::grid_provider_;
  using BaseType::scalar_ansatz_;
  using BaseType::scalar_test_;
  using BaseType::vector_ansatz_;
  using BaseType::vector_test_;
}; // struct ProductIntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class G>
using ProductIntegrandTest = Dune::GDT::Test::ProductIntegrandTest<G>;
TYPED_TEST_CASE(ProductIntegrandTest, Grids2D);

TYPED_TEST(ProductIntegrandTest, is_constructable)
{
  this->is_constructable();
}
TYPED_TEST(ProductIntegrandTest, evaluates_correctly_for_scalar_bases)
{
  this->evaluates_correctly_for_scalar_bases();
}

TYPED_TEST(ProductIntegrandTest, evaluates_correctly_for_vector_bases)
{
  this->evaluates_correctly_for_vector_bases();
}
