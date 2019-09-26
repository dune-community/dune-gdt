// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_INTEGRANDS_INTEGRANDS_HH
#define DUNE_GDT_TEST_INTEGRANDS_INTEGRANDS_HH

#include <array>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/type.hh>

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/fvector.hh>

#include <dune/xt/grid/gridprovider.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/xt/functions/generic/element-function.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct IntegrandTest : public ::testing::Test
{
  static_assert(XT::Grid::is_grid<G>::value, "");

  using GV = typename G::LeafGridView;
  using D = typename GV::ctype;
  static const constexpr size_t d = GV::dimension;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using LocalScalarBasisType = XT::Functions::GenericElementFunctionSet<E, 1, 1>;
  using DomainType = typename LocalScalarBasisType::DomainType;
  using ScalarRangeType = typename LocalScalarBasisType::RangeType;
  using ScalarJacobianType = typename LocalScalarBasisType::DerivativeRangeType;
  using LocalVectorBasisType = XT::Functions::GenericElementFunctionSet<E, 2, 1>;
  using VectorRangeType = typename LocalVectorBasisType::RangeType;
  using VectorJacobianType = typename LocalVectorBasisType::DerivativeRangeType;
  using MatrixType = typename XT::LA::Container<double, XT::LA::default_sparse_backend>::MatrixType;

  std::shared_ptr<XT::Grid::GridProvider<G>> grid_provider_;
  std::shared_ptr<LocalScalarBasisType> scalar_ansatz_;
  std::shared_ptr<LocalScalarBasisType> scalar_test_;
  std::shared_ptr<LocalVectorBasisType> vector_ansatz_;
  std::shared_ptr<LocalVectorBasisType> vector_test_;
  static constexpr bool is_simplex_grid_ = XT::Grid::is_uggrid<G>::value || XT::Grid::is_simplex_alugrid<G>::value;

  virtual std::shared_ptr<XT::Grid::GridProvider<G>> make_grid()
  {
    return std::make_shared<XT::Grid::GridProvider<G>>(
        XT::Grid::make_cube_grid<G>(XT::Common::from_string<FieldVector<double, d>>("[0 0 0 0]"),
                                    XT::Common::from_string<FieldVector<double, d>>("[3 1 1 1]"),
                                    XT::Common::from_string<std::array<unsigned int, d>>("[9 2 2 2]")));
  }

  void SetUp() override
  {
    grid_provider_ = make_grid();
    // {x, x^2 y}
    scalar_ansatz_ = std::make_shared<LocalScalarBasisType>(
        /*size = */ 2,
        /*ord = */ 3,
        /*evaluate = */
        [](const DomainType& x, std::vector<ScalarRangeType>& ret, const XT::Common::Parameter&) {
          ret = {{x[0]}, {std::pow(x[0], 2) * x[1]}};
        },
        /*param_type = */ XT::Common::ParameterType{},
        /*jacobian = */
        [](const DomainType& x, std::vector<ScalarJacobianType>& ret, const XT::Common::Parameter&) {
          // jacobians both only have a single row
          ret[0][0] = {1., 0.};
          ret[1][0] = {2. * x[0] * x[1], std::pow(x[0], 2)};
        });
    // {y, x y^3}
    scalar_test_ = std::make_shared<LocalScalarBasisType>(
        /*size = */ 2,
        /*ord = */ 4,
        /*evaluate = */
        [](const DomainType& x, std::vector<ScalarRangeType>& ret, const XT::Common::Parameter&) {
          ret = {{x[1]}, {x[0] * std::pow(x[1], 3)}};
        },
        /*param_type = */ XT::Common::ParameterType{},
        /*jacobian = */
        [](const DomainType& x, std::vector<ScalarJacobianType>& ret, const XT::Common::Parameter&) {
          // jacobians both only have a single row
          ret[0][0] = {0., 1.};
          ret[1][0] = {std::pow(x[1], 3), 3 * x[0] * std::pow(x[1], 2)};
        });
    // { (x,y)^T, (x^2, y^2)^T}
    vector_ansatz_ = std::make_shared<LocalVectorBasisType>(
        /*size = */ 2,
        /*ord = */ 2,
        /*evaluate = */
        [](const DomainType& x, std::vector<VectorRangeType>& ret, const XT::Common::Parameter&) {
          ret = {{x[0], x[1]}, {std::pow(x[0], 2), std::pow(x[1], 2)}};
        },
        /*param_type = */ XT::Common::ParameterType{},
        /*jacobian = */
        [](const DomainType& x, std::vector<VectorJacobianType>& ret, const XT::Common::Parameter&) {
          // jacobian of first function
          ret[0][0] = {1., 0.};
          ret[0][1] = {0., 1.};
          // jacobian of second function
          ret[1][0] = {2 * x[0], 0.};
          ret[1][1] = {0., 2 * x[1]};
        });
    // { (x,y)^T, (x^2, y^2)^T}
    vector_test_ = std::make_shared<LocalVectorBasisType>(
        /*size = */ 2,
        /*ord = */ 2,
        /*evaluate = */
        [](const DomainType& x, std::vector<VectorRangeType>& ret, const XT::Common::Parameter&) {
          ret = {{1, 2}, {x[0] * x[1], x[0] + x[1]}};
        },
        /*param_type = */ XT::Common::ParameterType{},
        /*jacobian = */
        [](const DomainType& x, std::vector<VectorJacobianType>& ret, const XT::Common::Parameter&) {
          // jacobian of first function
          ret[0][0] = {0., 0.};
          ret[0][1] = {0., 0.};
          // jacobian of second function
          ret[1][0] = {x[1], x[0]};
          ret[1][1] = {1., 1.};
        });
  } // ... SetUp(...)

  virtual void is_constructable() = 0;
}; // struct IntegrandTest


} // namespace Test
} // namespace GDT
} // namespace Dune


using Grids2D = ::testing::Types<YASP_2D_EQUIDISTANT_OFFSET
#if HAVE_DUNE_ALUGRID
                                 ,
                                 ALU_2D_SIMPLEX_CONFORMING,
                                 ALU_2D_SIMPLEX_NONCONFORMING,
                                 ALU_2D_CUBE
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
                                 ,
                                 UG_2D
#endif
#if HAVE_ALBERTA
                                 ,
                                 ALBERTA_2D
#endif
                                 >;

DUNE_XT_COMMON_TYPENAME(YASP_2D_EQUIDISTANT_OFFSET)
#if HAVE_DUNE_ALUGRID
DUNE_XT_COMMON_TYPENAME(ALU_2D_SIMPLEX_CONFORMING)
DUNE_XT_COMMON_TYPENAME(ALU_2D_SIMPLEX_NONCONFORMING)
DUNE_XT_COMMON_TYPENAME(ALU_2D_CUBE)
#endif
#if HAVE_DUNE_UGGRID || HAVE_UG
DUNE_XT_COMMON_TYPENAME(UG_2D)
#endif
#if HAVE_ALBERTA
DUNE_XT_COMMON_TYPENAME(ALBERTA_2D)
#endif


#endif // DUNE_GDT_TEST_INTEGRANDS_INTEGRANDS_HH
