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

#include <dune/xt/common/matrix.hh>

#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/local/integrands/symmetric_elliptic.hh>

#include <dune/gdt/test/integrands/integrands.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct SymmetricEllipticIntegrandTest : public IntegrandTest<G>
{
  using BaseType = IntegrandTest<G>;
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::E;
  using typename BaseType::GV;
  using typename BaseType::LocalVectorBasisType;
  using typename BaseType::VectorJacobianType;
  using typename BaseType::VectorRangeType;
  using VectorIntegrandType = LocalSymmetricEllipticIntegrand<E>;

  virtual void SetUp() override
  {
    BaseType::SetUp();
    diffusion_factor_ = std::make_shared<XT::Functions::GenericGridFunction<E, 1>>(
        2, [](const E&) {}, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    // { (x,x+y)^T, (x^2, x^2 + y^2)^T}
    vector_ansatz_ = std::make_shared<LocalVectorBasisType>(
        /*size = */ 2,
        /*ord = */ 2,
        /*evaluate = */
        [](const DomainType& x, std::vector<VectorRangeType>& ret, const XT::Common::Parameter&) {
          ret = {{x[0], x[0] + x[1]}, {std::pow(x[0], 2), std::pow(x[0], 2) + std::pow(x[1], 2)}};
        },
        /*param_type = */ XT::Common::ParameterType{},
        /*jacobian = */
        [](const DomainType& x, std::vector<VectorJacobianType>& ret, const XT::Common::Parameter&) {
          // jacobian of first function
          ret[0][0] = {1., 0.};
          ret[0][1] = {1., 1.};
          // jacobian of second function
          ret[1][0] = {2 * x[0], 0.};
          ret[1][1] = {2 * x[0], 2 * x[1]};
        });
  }

  virtual void is_constructable() override final
  {
    VectorIntegrandType vector_integrand1;
    VectorIntegrandType vector_integrand2(1.);
    const XT::Functions::GenericFunction<d, 1> scalar_function(
        2, [](const DomainType& x, const XT::Common::Parameter&) { return x[0] * x[1]; });
    VectorIntegrandType vector_integrand3(scalar_function);
    VectorIntegrandType vector_integrand4(*diffusion_factor_);
    DUNE_UNUSED_PARAMETER(vector_integrand1);
    DUNE_UNUSED_PARAMETER(vector_integrand2);
    DUNE_UNUSED_PARAMETER(vector_integrand3);
    DUNE_UNUSED_PARAMETER(vector_integrand4);
  }

  virtual void evaluates_correctly()
  {
    VectorIntegrandType integrand(*diffusion_factor_);
    const auto element = *(grid_provider_->leaf_view().template begin<0>());
    integrand.bind(element);
    const auto integrand_order = integrand.order(*vector_test_, *vector_ansatz_);
    EXPECT_EQ(4, integrand_order);
    DynamicMatrix<D> result(2, 2, 0.);
    for (const auto& quadrature_point : Dune::QuadratureRules<D, d>::rule(element.geometry().type(), integrand_order)) {
      const auto& x = quadrature_point.position();
      integrand.evaluate(*vector_test_, *vector_ansatz_, x, result);
      DynamicMatrix<D> expected_result{
          {0., 0.}, {x[1] + 0.5 * x[0] + 1.5, 2 * x[0] * x[1] + std::pow(x[0], 2) + x[0] + 2 * x[1]}};
      expected_result *= x[0] * x[1];
      for (size_t ii = 1; ii < 2; ++ii)
        for (size_t jj = 0; jj < 2; ++jj) {
          EXPECT_DOUBLE_EQ(expected_result[ii][jj], result[ii][jj]);
        }
    }
  }

  using BaseType::grid_provider_;
  using BaseType::scalar_ansatz_;
  using BaseType::scalar_test_;
  using BaseType::vector_ansatz_;
  using BaseType::vector_test_;
  std::shared_ptr<XT::Functions::GenericGridFunction<E, 1>> diffusion_factor_;
}; // struct SymmetricEllipticIntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


template <class G>
using SymmetricEllipticIntegrandTest = Dune::GDT::Test::SymmetricEllipticIntegrandTest<G>;
TYPED_TEST_CASE(SymmetricEllipticIntegrandTest, Grids2D);

TYPED_TEST(SymmetricEllipticIntegrandTest, is_constructable)
{
  this->is_constructable();
}

TYPED_TEST(SymmetricEllipticIntegrandTest, evaluates_correctly)
{
  this->evaluates_correctly();
}
