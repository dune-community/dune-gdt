// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#include <dune/xt/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/gdt/local/integrands/div.hh>

#include <dune/gdt/test/integrands/integrands.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct DivIntegrandTest : public IntegrandTest<G>
{
  using BaseType = IntegrandTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::MatrixType;
  using typename BaseType::VectorJacobianType;
  using TestDivIntegrandType = LocalElementAnsatzValueTestDivProductIntegrand<E>;
  using AnsatzDivIntegrandType = LocalElementAnsatzDivTestValueProductIntegrand<E>;

  virtual void is_constructable() override final
  {
    TestDivIntegrandType test_div_integrand1;
    TestDivIntegrandType test_div_integrand2(1.);
    const XT::Functions::GenericGridFunction<E, 1> scalar_grid_function(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    TestDivIntegrandType test_div_integrand3(scalar_grid_function);
    const XT::Functions::GenericFunction<d, 1> scalar_function(
        2, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    TestDivIntegrandType test_div_integrand4(scalar_function);
    DUNE_UNUSED_PARAMETER(test_div_integrand1);
    DUNE_UNUSED_PARAMETER(test_div_integrand2);
    DUNE_UNUSED_PARAMETER(test_div_integrand3);
    DUNE_UNUSED_PARAMETER(test_div_integrand4);
    AnsatzDivIntegrandType ansatz_div_integrand1;
    AnsatzDivIntegrandType ansatz_div_integrand2(1.);
    AnsatzDivIntegrandType ansatz_div_integrand3(scalar_grid_function);
    AnsatzDivIntegrandType ansatz_div_integrand4(scalar_function);
    DUNE_UNUSED_PARAMETER(ansatz_div_integrand1);
    DUNE_UNUSED_PARAMETER(ansatz_div_integrand2);
    DUNE_UNUSED_PARAMETER(ansatz_div_integrand3);
    DUNE_UNUSED_PARAMETER(ansatz_div_integrand4);
  }

  virtual void evaluates_correctly()
  {
    const XT::Functions::GenericGridFunction<E, 1> inducing_function(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    TestDivIntegrandType test_div_integrand(inducing_function);
    AnsatzDivIntegrandType ansatz_div_integrand(inducing_function);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    test_div_integrand.bind(element);
    ansatz_div_integrand.bind(element);
    const auto test_div_integrand_order = test_div_integrand.order(*vector_test_, *scalar_ansatz_);
    const auto ansatz_div_integrand_order = ansatz_div_integrand.order(*scalar_test_, *vector_ansatz_);
    DynamicMatrix<D> test_div_result(2, 2, 0.);
    DynamicMatrix<D> ansatz_div_result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(
             element.type(), std::max(test_div_integrand_order, ansatz_div_integrand_order))) {
      const auto& x = quadrature_point.position();
      test_div_integrand.evaluate(*vector_test_, *scalar_ansatz_, x, test_div_result);
      ansatz_div_integrand.evaluate(*scalar_test_, *vector_ansatz_, x, ansatz_div_result);
      DynamicMatrix<D> expected_result_test_div{{0, 0}, {x[0] * (x[1] + 1), std::pow(x[0], 2) * x[1] * (x[1] + 1)}};
      DynamicMatrix<D> expected_result_ansatz_div{
          {2 * x[1], 2 * (x[0] + x[1]) * x[1]},
          {2 * x[0] * std::pow(x[1], 3), 2 * (x[0] + x[1]) * x[0] * std::pow(x[1], 3)}};
      expected_result_test_div *= x[0] * x[1];
      expected_result_ansatz_div *= x[0] * x[1];
      for (size_t ii = 0; ii < 2; ++ii) {
        for (size_t jj = 0; jj < 2; ++jj) {
          EXPECT_DOUBLE_EQ(expected_result_test_div[ii][jj], test_div_result[ii][jj]);
          EXPECT_DOUBLE_EQ(expected_result_ansatz_div[ii][jj], ansatz_div_result[ii][jj]);
        } // jj
      } // ii
    } // quadrature points
  } // ... evaluates_correctly()

  virtual void is_integrated_correctly()
  {
    TestDivIntegrandType test_div_integrand(1.);
    AnsatzDivIntegrandType ansatz_div_integrand(1.);
    const auto& grid_view = grid_provider_->leaf_view();
    const auto scalar_space = make_continuous_lagrange_space<1>(grid_view, /*polorder=*/2);
    const auto vector_space = make_continuous_lagrange_space<d>(grid_view, /*polorder=*/2);
    const auto m = scalar_space.mapper().size();
    const auto n = vector_space.mapper().size();
    MatrixType test_div_mat(n, m, make_element_sparsity_pattern(vector_space, scalar_space, grid_view));
    MatrixType ansatz_div_mat(m, n, make_element_sparsity_pattern(scalar_space, vector_space, grid_view));
    MatrixOperator<MatrixType, GV, 1, 1, d, 1> test_div_op(grid_view, scalar_space, vector_space, test_div_mat);
    MatrixOperator<MatrixType, GV, d, 1, 1, 1> ansatz_div_op(grid_view, vector_space, scalar_space, ansatz_div_mat);
    test_div_op.append(LocalElementIntegralBilinearForm<E, d, 1, double, double, 1, 1>(test_div_integrand));
    test_div_op.assemble(true);
    ansatz_div_op.append(LocalElementIntegralBilinearForm<E, 1, 1, double, double, d, 1>(ansatz_div_integrand));
    ansatz_div_op.assemble(true);
    EXPECT_TRUE(XT::Common::FloatCmp::eq(test_div_mat, XT::Common::transposed(ansatz_div_mat), 1e-14, 1e-14));
    const auto mat_data_ptr = XT::Common::serialize_rowwise(test_div_mat);
    const auto min_entry = *std::min_element(mat_data_ptr.get(), mat_data_ptr.get() + n * m);
    const auto max_entry = *std::max_element(mat_data_ptr.get(), mat_data_ptr.get() + n * m);
    const auto square_sum = std::accumulate(
        mat_data_ptr.get(), mat_data_ptr.get() + n * m, 0., [](const auto& a, const auto& b) { return a + b * b; });
    EXPECT_NEAR((is_simplex_grid_ ? -0.133333333333333 : -0.177777777777778), min_entry, 1e-13);
    EXPECT_NEAR((is_simplex_grid_ ? 0.133333333333333 : 0.177777777777778), max_entry, 1e-13);
    EXPECT_NEAR((is_simplex_grid_ ? 4.515925925925922 : 4.589465020576143), square_sum, 1e-13);
    // std::cout << XT::Common::to_string(test_div_mat, 15) << std::endl;
  }

  using BaseType::grid_provider_;
  using BaseType::is_simplex_grid_;
  using BaseType::scalar_ansatz_;
  using BaseType::scalar_test_;
  using BaseType::vector_ansatz_;
  using BaseType::vector_test_;
}; // struct DivIntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class G>
using DivIntegrandTest = Dune::GDT::Test::DivIntegrandTest<G>;
TYPED_TEST_CASE(DivIntegrandTest, Grids2D);

TYPED_TEST(DivIntegrandTest, is_constructable)
{
  this->is_constructable();
}

TYPED_TEST(DivIntegrandTest, evaluates_correctly)
{
  this->evaluates_correctly();
}

TYPED_TEST(DivIntegrandTest, integrates_correctly)
{
  this->is_integrated_correctly();
}
