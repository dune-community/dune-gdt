// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH
#define DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH

#include <dune/xt/common/test/common.hh>
#include <dune/xt/common/test/gtest/gtest.h>

#include <dune/xt/la/container/istl.hh>

#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>

#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/tools/dirichlet-constraints.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class GV>
struct StokesDirichletProblem
{
  static_assert(XT::Grid::is_view<GV>::value, "");
  static_assert(GV::dimension == 2, "");

  static const constexpr size_t d = GV::dimension;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using G = typename GV::Grid;
  using RangeField = double;

  StokesDirichletProblem(
      std::shared_ptr<const XT::Functions::GridFunctionInterface<E, 1, 1>> diffusion_factor_in =
          std::make_shared<const XT::Functions::ConstantGridFunction<E, 1, 1>>(1., "trivial diffusion factor"),
      std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, 1>> rhs_in =
          std::make_shared<const XT::Functions::ConstantGridFunction<E, d, 1>>(0., "zero rhs"),
      std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, 1>> dirichlet_in =
          std::make_shared<const XT::Functions::ConstantGridFunction<E, d, 1>>(0., "dirichlet zero boundary values"))
    : diffusion_factor_(diffusion_factor_in)
    , diffusion_tensor_(std::make_shared<const XT::Functions::ConstantGridFunction<E, d, d>>(
          XT::LA::eye_matrix<FieldMatrix<double, d, d>>(d, d)))
    , rhs_(rhs_in)
    , dirichlet_(dirichlet_in)
    , boundary_info_()
    , grid_(XT::Grid::make_cube_grid<G>(-1, 1, 10))
    , grid_view_(grid_.leaf_view())
  {}

  const GV& grid_view()
  {
    return grid_view_;
  } // ... make_initial_grid(...)

  const XT::Functions::GridFunctionInterface<E, 1, 1>& diffusion_factor()
  {
    return *diffusion_factor_;
  } // ... make_initial_grid(...)

  const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor()
  {
    return *diffusion_tensor_;
  } // ... make_initial_grid(...)

  const XT::Functions::GridFunctionInterface<E, d, 1>& rhs()
  {
    return *rhs_;
  } // ... make_initial_grid(...)

  const XT::Functions::GridFunctionInterface<E, d, 1>& dirichlet()
  {
    return *dirichlet_;
  } // ... make_initial_grid(...)

  const XT::Grid::BoundaryInfo<I>& boundary_info()
  {
    return boundary_info_;
  } // ... make_initial_grid(...)

  std::shared_ptr<const XT::Functions::GridFunctionInterface<E, 1, 1, RangeField>> diffusion_factor_;
  std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, d, RangeField>> diffusion_tensor_;
  std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, 1, RangeField>> rhs_;
  std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, 1, RangeField>> dirichlet_;
  const XT::Grid::AllDirichletBoundaryInfo<I> boundary_info_;
  const XT::Grid::GridProvider<G> grid_;
  const GV grid_view_;
}; // class StokesDirichletProblem


template <class G>
class StokesDirichletTest : public ::testing::Test
{
  using GV = typename G::LeafGridView;
  static const size_t d = GV::dimension;

  using Matrix = XT::LA::IstlRowMajorSparseMatrix<double>;
  using Vector = XT::LA::IstlDenseVector<double>;
  using VelocitySpace = ContinuousLagrangeSpace<GV, d>;
  using PressureSpace = ContinuousLagrangeSpace<GV, 1>;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using RangeField = double;

public:
  StokesDirichletTest(StokesDirichletProblem<GV> problem)
    : problem_(problem)
  {}

  void run()
  {
    const auto& grid_view = problem_.grid_view();
    const VelocitySpace velocity_space(grid_view, 2);
    const PressureSpace pressure_space(grid_view, 1);
    const size_t size_u = velocity_space.mapper().size();
    const size_t size_p = pressure_space.mapper().size();
    Matrix A(size_u, size_u, make_element_sparsity_pattern(velocity_space, velocity_space, grid_view));
    Matrix B(size_u, size_p, make_element_sparsity_pattern(velocity_space, pressure_space, grid_view));
    MatrixOperator<Matrix, GV, d> A_operator(grid_view, velocity_space, velocity_space, A);
    MatrixOperator<Matrix, GV, 1, 1, d> B_operator(grid_view, pressure_space, velocity_space, B);
    A_operator.append(LocalElementIntegralBilinearForm<E, d>(
        LocalEllipticIntegrand<E, d>(problem_.diffusion_factor(), problem_.diffusion_tensor())));
    B_operator.append(
        LocalElementIntegralBilinearForm<E, d, 1, RangeField, RangeField, 1>(LocalElementProductTestDivIntegrand<E>()));
    auto rhs_u = make_vector_functional<Vector>(velocity_space);
    rhs_u.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), problem_.rhs())));
    A_operator.append(rhs_u);
    // Dirichlet constrainst for u
    DirichletConstraints<I, VelocitySpace> dirichlet_constraints(problem_.boundary_info(), velocity_space);
    A_operator.append(dirichlet_constraints);
    // assemble everything in one grid walk
    problem_.diffusion_factor().visualize(grid_view, "testvis");
    A_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    B_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    dirichlet_constraints.apply(A, rhs_u.vector());
  }

  StokesDirichletProblem<GV> problem_;
}; // class StrokesDirichletTest

template <class G>
class StokesTestcase1 : public StokesDirichletTest<G>
{
  using BaseType = StokesDirichletTest<G>;

public:
  StokesTestcase1()
    : BaseType(StokesDirichletProblem<typename G::LeafGridView>())
  {}
}; // struct StokesTestcase1


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH
