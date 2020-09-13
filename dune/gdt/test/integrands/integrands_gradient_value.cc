// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/geometry/quadraturerules.hh>

#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/local/integrands/gradient-value.hh>

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
  using typename BaseType::MatrixType;
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
    [[maybe_unused]] ScalarIntegrandType scalar_integrand1(vector_grid_function);
    [[maybe_unused]] ScalarIntegrandType scalar_integrand2(vector_function);
    [[maybe_unused]] VectorIntegrandType vector_integrand1(vector_grid_function);
    [[maybe_unused]] VectorIntegrandType vector_integrand2(vector_function);
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
    DynamicMatrix<D> result(2, 2, 0.);
    DynamicMatrix<D> result2(2, 2, 0.);
    for (const auto& quadrature_point :
         Dune::QuadratureRules<D, d>::rule(element.type(), std::max(integrand_order, integrand_order2))) {
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
    DynamicMatrix<D> result(2, 2, 0.);
    DynamicMatrix<D> result2(2, 2, 0.);
    for (const auto& quadrature_point :
         Dune::QuadratureRules<D, d>::rule(element.type(), std::max(integrand_order, integrand_order2))) {
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

  virtual void is_integrated_correctly()
  {
    const XT::Functions::GenericFunction<d, d> vector_function(
        0, [](const DomainType& /*x*/, const XT::Common::Parameter&) {
          return FieldVector<D, d>{{1., 1.}};
        });
    VectorIntegrandType ansatz_grad_integrand(vector_function);
    VectorIntegrandTestGradType test_grad_integrand(vector_function);
    const auto& grid_view = grid_provider_->leaf_view();
    const auto space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    const auto n = space.mapper().size();
    MatrixType test_grad_mat(n, n, make_element_sparsity_pattern(space, space, grid_view));
    MatrixType ansatz_grad_mat(n, n, make_element_sparsity_pattern(space, space, grid_view));
    MatrixOperator<MatrixType, GV, d, 1, d, 1> test_grad_op(grid_view, space, space, test_grad_mat);
    MatrixOperator<MatrixType, GV, d, 1, d, 1> ansatz_grad_op(grid_view, space, space, ansatz_grad_mat);
    test_grad_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, d, 1>{test_grad_integrand});
    test_grad_op.assemble(true);
    ansatz_grad_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, d, 1>{ansatz_grad_integrand});
    ansatz_grad_op.assemble(true);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(test_grad_mat, XT::Common::transposed(ansatz_grad_mat), 1e-14, 1e-14));
    const auto mat_data_ptr = XT::Common::serialize_rowwise(test_grad_mat);
    const auto min_entry = *std::min_element(mat_data_ptr.get(), mat_data_ptr.get() + n * n);
    const auto max_entry = *std::max_element(mat_data_ptr.get(), mat_data_ptr.get() + n * n);
    const auto square_sum = std::accumulate(
        mat_data_ptr.get(), mat_data_ptr.get() + n * n, 0., [](const auto& a, const auto& b) { return a + b * b; });
    EXPECT_NEAR((is_simplex_grid_ ? -0.133333333333333 : -0.177777777777778), min_entry, 1e-13);
    EXPECT_NEAR((is_simplex_grid_ ? 0.133333333333333 : 0.177777777777778), max_entry, 1e-13);
    EXPECT_NEAR((is_simplex_grid_ ? 5.208148148148139 : 9.178930041152277), square_sum, 1e-13);
    // std::cout << XT::Common::to_string(test_grad_mat, 15) << std::endl;
  }

  using BaseType::grid_provider_;
  using BaseType::is_simplex_grid_;
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
TYPED_TEST_SUITE(GradientValueIntegrandTest, Grids2D);

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

TYPED_TEST(GradientValueIntegrandTest, is_integrated_correctly)
{
  this->is_integrated_correctly();
}
