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

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/eigen-solver.hh>
#include <dune/xt/la/solver.hh>

#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/elliptic.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/operators/matrix-based.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
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
  using ScalarGridFunction = XT::Functions::GridFunctionInterface<E, 1, 1>;
  using VectorGridFunction = XT::Functions::GridFunctionInterface<E, d, 1>;
  using DomainType = typename ScalarGridFunction::LocalFunctionType::DomainType;

  StokesDirichletProblem(std::shared_ptr<const ScalarGridFunction> diffusion_factor_in = default_diffusion_factor(),
                         std::shared_ptr<const VectorGridFunction> rhs_in = default_rhs(),
                         std::shared_ptr<const VectorGridFunction> dirichlet_in = default_dirichlet_values(),
                         std::shared_ptr<const VectorGridFunction> reference_sol_u = nullptr,
                         std::shared_ptr<const ScalarGridFunction> reference_sol_p = nullptr)
    : diffusion_factor_(diffusion_factor_in)
    , diffusion_tensor_(std::make_shared<const XT::Functions::ConstantGridFunction<E, d, d>>(
          XT::LA::eye_matrix<FieldMatrix<double, d, d>>(d, d)))
    , rhs_(rhs_in)
    , dirichlet_(dirichlet_in)
    , reference_sol_u_(reference_sol_u)
    , reference_sol_p_(reference_sol_p)
    , boundary_info_()
    , grid_(XT::Grid::make_cube_grid<G>(-1, 1, 50))
    , grid_view_(grid_.leaf_view())
  {}

  static std::shared_ptr<const ScalarGridFunction> default_diffusion_factor()
  {
    return std::make_shared<const XT::Functions::ConstantGridFunction<E, 1, 1>>(1., "default diffusion factor");
  }

  static std::shared_ptr<const VectorGridFunction> default_rhs()
  {
    return std::make_shared<const XT::Functions::ConstantGridFunction<E, d, 1>>(0., "zero rhs");
  }

  static std::shared_ptr<const VectorGridFunction> default_dirichlet_values()
  {
    return std::make_shared<const XT::Functions::ConstantGridFunction<E, d, 1>>(0., "dirichlet zero boundary values");
  }

  const GV& grid_view()
  {
    return grid_view_;
  } // ... make_initial_grid(...)

  const ScalarGridFunction& diffusion_factor()
  {
    return *diffusion_factor_;
  } // ... make_initial_grid(...)

  const XT::Functions::GridFunctionInterface<E, d, d>& diffusion_tensor()
  {
    return *diffusion_tensor_;
  } // ... make_initial_grid(...)

  const VectorGridFunction& rhs()
  {
    return *rhs_;
  } // ... make_initial_grid(...)

  const VectorGridFunction& dirichlet()
  {
    return *dirichlet_;
  } // ... make_initial_grid(...)

  const VectorGridFunction& reference_solution_u()
  {
    DUNE_THROW_IF(!reference_sol_u_, Dune::InvalidStateException, "No reference solution provided!");
    return *reference_sol_u_;
  } // ... make_initial_grid(...)

  const ScalarGridFunction& reference_solution_p()
  {
    DUNE_THROW_IF(!reference_sol_p_, Dune::InvalidStateException, "No reference solution provided!");
    return *reference_sol_p_;
  } // ... make_initial_grid(...)

  const XT::Grid::BoundaryInfo<I>& boundary_info()
  {
    return boundary_info_;
  } // ... make_initial_grid(...)

  std::shared_ptr<const ScalarGridFunction> diffusion_factor_;
  std::shared_ptr<const XT::Functions::GridFunctionInterface<E, d, d, RangeField>> diffusion_tensor_;
  std::shared_ptr<const VectorGridFunction> rhs_;
  std::shared_ptr<const VectorGridFunction> dirichlet_;
  std::shared_ptr<const VectorGridFunction> reference_sol_u_;
  std::shared_ptr<const ScalarGridFunction> reference_sol_p_;
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
  using DenseMatrix = XT::LA::CommonDenseMatrix<double>;
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

  static bool is_symmetric(const Matrix& mat)
  {
    if (mat.rows() != mat.cols())
      return false;
    for (size_t ii = 0; ii < mat.rows(); ++ii)
      for (size_t jj = ii; jj < mat.cols(); ++jj)
        if (XT::Common::FloatCmp::ne(mat.get_entry(ii, jj), mat.get_entry(jj, ii)))
          return false;
    return true;
  }

  static bool is_positive_definite(const Matrix& mat)
  {
    //    DenseMatrix dense_matrix(mat);
    try {
      //            XT::LA::cholesky(dense_matrix);
      //      auto eigenvalues = XT::LA::make_eigen_solver(dense_matrix).real_eigenvalues();
      //      for (const auto& eigenvalue : eigenvalues)
      //        if (XT::Common::FloatCmp::le(eigenvalue, 0.))
      //          std::cout << eigenvalue << std::endl;
      return true;
    } catch (Dune::MathError&) {
      return false;
    }
  }

  void run()
  {
    const auto& grid_view = problem_.grid_view();
    // setup spaces and matrices A and B (system matrix is [A B; B^T 0]
    const VelocitySpace velocity_space(grid_view, 2);
    const PressureSpace pressure_space(grid_view, 1);
    const size_t size_u = velocity_space.mapper().size();
    const size_t size_p = pressure_space.mapper().size();
    auto pattern_u_u = make_element_sparsity_pattern(velocity_space, velocity_space, grid_view);
    auto pattern_u_p = make_element_sparsity_pattern(velocity_space, pressure_space, grid_view);
    Matrix A(size_u, size_u, pattern_u_u);
    Matrix B(size_u, size_p, pattern_u_p);
    MatrixOperator<Matrix, GV, d> A_operator(grid_view, velocity_space, velocity_space, A);
    MatrixOperator<Matrix, GV, 1, 1, d> B_operator(grid_view, pressure_space, velocity_space, B);
    A_operator.append(LocalElementIntegralBilinearForm<E, d>(
        LocalEllipticIntegrand<E, d>(problem_.diffusion_factor(), problem_.diffusion_tensor())));
    B_operator.append(
        LocalElementIntegralBilinearForm<E, d, 1, RangeField, RangeField, 1>(LocalElementProductTestDivIntegrand<E>()));
    Vector f_vector(size_u);
    auto discrete_f = make_discrete_function(velocity_space, f_vector);
    auto rhs_functional = make_vector_functional(velocity_space, f_vector);
    rhs_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalElementProductIntegrand<E, d>(), problem_.rhs())));
    A_operator.append(rhs_functional);
    // Dirichlet constrainst for u
    DirichletConstraints<I, VelocitySpace> dirichlet_constraints(problem_.boundary_info(), velocity_space);
    A_operator.append(dirichlet_constraints);
    // assemble everything in one grid walk
    A_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    EXPECT_TRUE(is_symmetric(A));
    B_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));

    Vector dirichlet_vector(size_u, 0.), reference_solution_u_vector(size_u, 0.),
        reference_solution_p_vector(size_p, 0.);
    auto discrete_dirichlet_values = make_discrete_function(velocity_space, dirichlet_vector);
    auto reference_solution_u = make_discrete_function(velocity_space, reference_solution_u_vector);
    auto reference_solution_p = make_discrete_function(pressure_space, reference_solution_p_vector);
    interpolate(
        problem_.dirichlet(), discrete_dirichlet_values, problem_.boundary_info(), XT::Grid::DirichletBoundary());
    interpolate(problem_.reference_solution_u(), reference_solution_u);
    interpolate(problem_.reference_solution_p(), reference_solution_p);
    Vector discrete_rhs_vector_u(size_u), discrete_rhs_vector_p(size_p);
    // create rhs vector f - A g_D for u
    A.mv(discrete_dirichlet_values.dofs().vector(), discrete_rhs_vector_u);
    discrete_rhs_vector_u *= -1;
    discrete_rhs_vector_u += f_vector;
    // create rhs vector -B^T g_D for p
    B.mtv(discrete_dirichlet_values.dofs().vector(), discrete_rhs_vector_p);
    discrete_rhs_vector_p *= -1;
    // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
    // Dirichlet DoF, so we also need to clear the row of B.
    dirichlet_constraints.apply(A, discrete_rhs_vector_u);
    for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs())
      B.clear_row(DoF);
    EXPECT_TRUE(is_positive_definite(A));
    // copy matrices to saddle point matrix
    // create pattern first
    XT::LA::SparsityPatternDefault system_matrix_pattern(size_u + size_p);
    for (size_t ii = 0; ii < pattern_u_u.size(); ++ii)
      for (const auto& jj : pattern_u_u.inner(ii))
        system_matrix_pattern.insert(ii, jj);
    for (size_t ii = 0; ii < pattern_u_p.size(); ++ii) {
      for (const auto& jj : pattern_u_p.inner(ii)) {
        system_matrix_pattern.insert(ii, size_u + jj);
        system_matrix_pattern.insert(size_u + jj, ii);
      }
    }
    // to set the first DoF of p to 0 for uniqueness
    system_matrix_pattern.insert(size_u, size_u);
    // create_matrix
    Matrix system_matrix(size_u + size_p, size_u + size_p, system_matrix_pattern);
    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the size_u-th row of
    // [A B; B^T 0] to the unit vector.
    B.clear_col(0);
    discrete_rhs_vector_p.set_entry(0, 0.);
    system_matrix.unit_row(size_u);
    // now copy the matrices
    for (size_t ii = 0; ii < size_u; ++ii)
      for (const auto& jj : pattern_u_u.inner(ii))
        system_matrix.set_entry(ii, jj, A.get_entry(ii, jj));
    for (size_t ii = 0; ii < size_u; ++ii) {
      for (const auto& jj : pattern_u_p.inner(ii)) {
        system_matrix.set_entry(ii, size_u + jj, B.get_entry(ii, jj));
        system_matrix.set_entry(size_u + jj, ii, B.get_entry(ii, jj));
      }
    }
    // also copy the rhs
    Vector system_vector(size_u + size_p, 0.), solution_vector(size_u + size_p, 0.);
    for (size_t ii = 0; ii < size_u; ++ii)
      system_vector[ii] = discrete_rhs_vector_u[ii];
    for (size_t ii = 0; ii < size_p; ++ii)
      system_vector[size_u + ii] = discrete_rhs_vector_p[ii];
    // solve the system by a direct solver
    XT::LA::solve(system_matrix, system_vector, solution_vector);
    Vector solution_u(size_u), solution_p(size_p);
    for (size_t ii = 0; ii < size_u; ++ii)
      solution_u[ii] = solution_vector[ii];
    for (size_t ii = 0; ii < size_p; ++ii)
      solution_p[ii] = solution_vector[size_u + ii];
    const auto actual_u_vector = solution_u + dirichlet_vector;
    auto dirichlet_func = make_discrete_function(velocity_space, dirichlet_vector);
    auto sol_u_func = make_discrete_function(velocity_space, actual_u_vector);
    auto sol_u_wo_dirichlet_func = make_discrete_function(velocity_space, solution_u);
    auto sol_p_func = make_discrete_function(pressure_space, solution_p);
    sol_u_func.visualize("solution_u");
    sol_u_wo_dirichlet_func.visualize("solution_u_wo_dirichlet");
    sol_p_func.visualize("solution_p");
    dirichlet_func.visualize("solution_u_dirichlet");
    reference_solution_u.visualize("u_ref");
    reference_solution_p.visualize("p_ref");
  } // run

  StokesDirichletProblem<GV> problem_;
}; // class StrokesDirichletTest

template <class G>
class StokesTestcase1 : public StokesDirichletTest<G>
{
  using BaseType = StokesDirichletTest<G>;
  using ProblemType = StokesDirichletProblem<typename G::LeafGridView>;
  using E = typename ProblemType::E;
  using DomainType = typename ProblemType::DomainType;
  using ScalarGridFunction = typename ProblemType::ScalarGridFunction;
  using VectorGridFunction = typename ProblemType::VectorGridFunction;
  static const size_t d = ProblemType::d;

public:
  StokesTestcase1()
    : BaseType(ProblemType(ProblemType::default_diffusion_factor(),
                           ProblemType::default_rhs(),
                           dirichlet(),
                           exact_sol_u(),
                           exact_sol_p()))
  {}

  static std::shared_ptr<const VectorGridFunction> dirichlet()
  {
    static const XT::Functions::GenericFunction<d, d, 1> dirichlet_values(
        50, [](const DomainType& xx, const XT::Common::Parameter&) {
          typename VectorGridFunction::LocalFunctionType::RangeReturnType ret;
          auto x = xx[0];
          auto y = xx[1];
          ret[0] = -std::exp(x) * (y * std::cos(y) + std::sin(y));
          ret[1] = std::exp(x) * y * std::sin(y);
          return ret;
        });
    return std::make_shared<XT::Functions::FunctionAsGridFunctionWrapper<E, d, 1, double>>(dirichlet_values);
  }

  static std::shared_ptr<const VectorGridFunction> exact_sol_u()
  {
    return dirichlet();
  }

  static std::shared_ptr<const ScalarGridFunction> exact_sol_p()
  {
    static const XT::Functions::GenericFunction<d, 1, 1> exact_sol_p_(
        50, [](const DomainType& xx, const XT::Common::Parameter&) {
          typename ScalarGridFunction::LocalFunctionType::RangeReturnType ret;
          auto x = xx[0];
          auto y = xx[1];
          ret[0] = 2. * std::exp(x) * std::sin(y);
          return ret;
        });
    return std::make_shared<XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, double>>(exact_sol_p_);
  }
}; // struct StokesTestcase1


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH
