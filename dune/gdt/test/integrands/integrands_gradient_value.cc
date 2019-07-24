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

#include <dune/gdt/local/integrands/gradient_value.hh>

#include <dune/gdt/test/integrands/integrands.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct GradientValueIntegrandTest : public IntegrandTest<G>
{
  using BaseType = IntegrandTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::VectorJacobianType;
  using ScalarIntegrandType = LocalElementGradientValueIntegrand<E, 1>;
  using ScalarIntegrandTestGradType = LocalElementGradientValueIntegrand<E, 1, 1, D, D, D, true>;
  using VectorIntegrandType = LocalElementGradientValueIntegrand<E, d>;
  using VectorIntegrandTestGradType = LocalElementGradientValueIntegrand<E, d, 1, D, D, D, true>;

  virtual void is_constructable() override final
  {
    const XT::Functions::GenericGridFunction<E, d> vector_grid_function(
        2,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return FieldVector<D, d>{{x[0], x[0] * x[1]}};
        });
    const XT::Functions::GenericFunction<d, d> vector_function(2,
                                                               [](const DomainType& x, const XT::Common::Parameter&) {
                                                                 return FieldVector<D, d>{{x[0], x[0] * x[1]}};
                                                               });
    ScalarIntegrandType scalar_integrand1(vector_grid_function);
    ScalarIntegrandType scalar_integrand2(vector_function);
    DUNE_UNUSED_PARAMETER(scalar_integrand1);
    DUNE_UNUSED_PARAMETER(scalar_integrand2);
    VectorIntegrandType vector_integrand1(vector_grid_function);
    VectorIntegrandType vector_integrand2(vector_function);
    DUNE_UNUSED_PARAMETER(vector_integrand1);
    DUNE_UNUSED_PARAMETER(vector_integrand2);
  }

  virtual void evaluates_correctly_for_scalar_bases()
  {
    const XT::Functions::GenericGridFunction<E, d> vector_grid_function(
        2,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return FieldVector<D, d>{{x[0], x[0] * x[1]}};
        });
    ScalarIntegrandType scalar_integrand(vector_grid_function);
    ScalarIntegrandTestGradType scalar_integrand2(vector_grid_function);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    scalar_integrand.bind(element);
    scalar_integrand2.bind(element);
    const auto integrand_order = scalar_integrand.order(*scalar_test_, *scalar_ansatz_);
    const auto integrand_order2 = scalar_integrand2.order(*scalar_test_, *scalar_ansatz_);
    EXPECT_EQ(8, integrand_order);
    EXPECT_EQ(integrand_order, integrand_order2);
    DynamicMatrix<D> result(2, 2, 0.);
    DynamicMatrix<D> result2(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      scalar_integrand.evaluate(*scalar_test_, *scalar_ansatz_, x, result);
      scalar_integrand2.evaluate(*scalar_test_, *scalar_ansatz_, x, result2);
      DynamicMatrix<D> expected_result{
          {x[0] * x[1], (2 + x[0]) * std::pow(x[0] * x[1], 2)},
          {std::pow(x[0] * x[1], 2) * x[1], (2 * x[1] + x[0] * x[1]) * std::pow(x[0] * x[1], 3)}};
      DynamicMatrix<D> expected_result2{{std::pow(x[0], 2) * x[1], std::pow(x[0] * x[1], 2) * x[0]},
                                        {std::pow(x[1], 3) * (std::pow(x[0], 2) + 3 * std::pow(x[0], 3)),
                                         std::pow(x[1], 4) * (std::pow(x[0], 3) + 3 * std::pow(x[0], 4))}};
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < 2; ++jj) {
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
          EXPECT_DOUBLE_EQ(expected_result2[ii][jj], result2[ii][jj]);
        } // jj
      } // ii
    } // quad_points
  }

  virtual void evaluates_correctly_for_vector_bases()
  {
    const XT::Functions::GenericGridFunction<E, d> vector_grid_function(
        2,
        [](const E&) {},
        [](const DomainType& x, const XT::Common::Parameter&) {
          return FieldVector<D, d>{{x[0], x[0] * x[1]}};
        });
    VectorIntegrandType integrand(vector_grid_function);
    VectorIntegrandTestGradType integrand2(vector_grid_function);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    integrand.bind(element);
    integrand2.bind(element);
    const auto integrand_order = integrand.order(*vector_test_, *vector_ansatz_);
    const auto integrand_order2 = integrand2.order(*vector_test_, *vector_ansatz_);
    EXPECT_EQ(5, integrand_order);
    EXPECT_EQ(integrand_order, integrand_order2);
    DynamicMatrix<D> result(2, 2, 0.);
    DynamicMatrix<D> result2(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      integrand.evaluate(*vector_test_, *vector_ansatz_, x, result);
      integrand2.evaluate(*vector_test_, *vector_ansatz_, x, result2);
      DynamicMatrix<D> expected_result{{x[0] * (1 + 2 * x[1]), 2 * std::pow(x[0], 2) + 4 * x[0] * std::pow(x[1], 2)},
                                       {2 * std::pow(x[0], 2) * x[1] + x[0] * std::pow(x[1], 2),
                                        2 * std::pow(x[0], 3) * x[1] + 2 * x[0] * std::pow(x[1], 2) * (x[0] + x[1])}};
      DynamicMatrix<D> expected_result2{
          {0, 0},
          {x[1] * (x[0] + std::pow(x[0], 2) + std::pow(x[0], 3)) + x[0] * std::pow(x[1], 2),
           x[1] * (std::pow(x[0], 3) + std::pow(x[0], 4)) + x[0] * (std::pow(x[1], 2) + std::pow(x[1], 3))}};
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < 2; ++jj) {
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
          EXPECT_DOUBLE_EQ(expected_result2[ii][jj], result2[ii][jj]);
        } // jj
      } // ii
    } // quad_points
  }

  using BaseType::grid_provider_;
  using BaseType::scalar_ansatz_;
  using BaseType::scalar_test_;
  using BaseType::vector_ansatz_;
  using BaseType::vector_test_;
}; // struct GradientValueIntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class G>
using GradientValueIntegrandTest = Dune::GDT::Test::GradientValueIntegrandTest<G>;
TYPED_TEST_CASE(GradientValueIntegrandTest, Grids2D);

TYPED_TEST(GradientValueIntegrandTest, is_constructable)
{
  this->is_constructable();
}
TYPED_TEST(GradientValueIntegrandTest, evaluates_correctly_for_scalar_bases)
{
  this->evaluates_correctly_for_scalar_bases();
}

TYPED_TEST(GradientValueIntegrandTest, evaluates_correctly_for_vector_bases)
{
  this->evaluates_correctly_for_vector_bases();
}
