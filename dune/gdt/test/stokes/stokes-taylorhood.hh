// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH
#define DUNE_GDT_TEST_STOKES_STOKES_TAYLORHOOD_HH

#include <dune/xt/test/common.hh>
#include <dune/xt/test/common/float_cmp.hh>
#include <dune/xt/test/gtest/gtest.h>

#include <dune/xt/la/algorithms/cholesky.hh>
#include <dune/xt/la/container/istl.hh>
#include <dune/xt/la/solver/istl/saddlepoint.hh>

#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/xt/functions/base/function-as-grid-function.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/generic/grid-function.hh>

#include <dune/gdt/functionals/vector-based.hh>
#include <dune/gdt/interpolations/boundary.hh>
#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/local/bilinear-forms/integrals.hh>
#include <dune/gdt/local/integrands/conversion.hh>
#include <dune/gdt/local/integrands/div.hh>
#include <dune/gdt/local/integrands/laplace.hh>
#include <dune/gdt/local/integrands/product.hh>
#include <dune/gdt/local/functionals/integrals.hh>
#include <dune/gdt/norms.hh>
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
  using DiffusionTensor = XT::Functions::GridFunctionInterface<E, d, d, RangeField>;

  StokesDirichletProblem(std::shared_ptr<const DiffusionTensor> diffusion_in = default_diffusion,
                         std::shared_ptr<const VectorGridFunction> rhs_f_in = default_rhs_f(),
                         std::shared_ptr<const VectorGridFunction> rhs_g_in = default_rhs_g(),
                         std::shared_ptr<const VectorGridFunction> dirichlet_in = default_dirichlet_values(),
                         std::shared_ptr<const VectorGridFunction> reference_sol_u = nullptr,
                         std::shared_ptr<const ScalarGridFunction> reference_sol_p = nullptr)
    : diffusion_(diffusion_in)
    , rhs_f_(rhs_f_in)
    , rhs_g_(rhs_g_in)
    , dirichlet_(dirichlet_in)
    , reference_sol_u_(reference_sol_u)
    , reference_sol_p_(reference_sol_p)
    , boundary_info_()
    , grid_(XT::Grid::make_cube_grid<G>(-1, 1, 20))
    , grid_view_(grid_.leaf_view())
  {}

  static std::shared_ptr<const DiffusionTensor> default_diffusion()
  {
    return std::make_shared<const XT::Functions::ConstantGridFunction<E, d, d>>(
        XT::LA::eye_matrix<FieldMatrix<double, d, d>>(d, d), "isotropic_diffusion");
  }

  static std::shared_ptr<const VectorGridFunction> default_rhs_f()
  {
    return std::make_shared<const XT::Functions::ConstantGridFunction<E, d, 1>>(0., "zero rhs");
  }

  static std::shared_ptr<const VectorGridFunction> default_rhs_g()
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

  const XT::Functions::GridFunctionInterface<E, d, d>& diffusion()
  {
    return *diffusion_;
  } // ... make_initial_grid(...)

  const VectorGridFunction& rhs_f()
  {
    return *rhs_f_;
  } // ... make_initial_grid(...)

  const VectorGridFunction& rhs_g()
  {
    return *rhs_g_;
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

  std::shared_ptr<const DiffusionTensor> diffusion_;
  std::shared_ptr<const VectorGridFunction> rhs_f_;
  std::shared_ptr<const VectorGridFunction> rhs_g_;
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
        if (XT::Common::FloatCmp::ne(mat.get_entry(ii, jj), mat.get_entry(jj, ii))) {
          std::cerr << "mat.rows() = " << mat.rows() << ", mat.cols() = " << mat.cols() << std::endl;
          std::cerr << "not symmetric for indices " << ii << ", " << jj << " with values " << mat.get_entry(ii, jj)
                    << ", " << mat.get_entry(jj, ii) << std::endl;
          return false;
        }
    return true;
  }

  static bool is_positive_definite(const Matrix& mat)
  {
    DenseMatrix dense_matrix(mat);
    try {
      XT::LA::cholesky(dense_matrix);
      return true;
    } catch (Dune::MathError&) {
      return false;
    }
  }

  void run(const int velocity_order, const double expected_error_u, const double expected_error_p)
  {
    const auto& grid_view = problem_.grid_view();
    // Setup spaces and matrices and vectors
    // Equations are
    // \int \nabla u \nabla v - \int p div v = \int ff v
    // \int (div u) q = \int gg q
    // System is [A B; B^T C] [u; p] = [f; g]
    // Dimensions are: A: n x n, B: n x m, C: m x m, u: n, f: n, p: m, g: m
    const VelocitySpace velocity_space(grid_view, velocity_order);
    const PressureSpace pressure_space(grid_view, velocity_order - 1);
    const size_t m = velocity_space.mapper().size();
    const size_t n = pressure_space.mapper().size();
    auto pattern_A = make_element_sparsity_pattern(velocity_space, velocity_space, grid_view);
    auto pattern_B = make_element_sparsity_pattern(velocity_space, pressure_space, grid_view);
    auto pattern_C = make_element_sparsity_pattern(pressure_space, pressure_space, grid_view);
    Matrix A(m, m, pattern_A);
    Matrix B(m, n, pattern_B);
    MatrixOperator<Matrix, GV, d> A_operator(grid_view, velocity_space, velocity_space, A);
    MatrixOperator<Matrix, GV, 1, 1, d> B_operator(grid_view, pressure_space, velocity_space, B);
    // calculate A_{ij} as \int \nabla v_i \nabla v_j
    A_operator.append(LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(problem_.diffusion())));
    // calculate B_{ij} as \int -\nabla p_i div(v_j)
    B_operator.append(LocalElementIntegralBilinearForm<E, d, 1, RangeField, RangeField, 1>(
        LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.)));
    // calculate rhs f as \int ff v and the integrated pressure space basis \int q_i
    Vector f_vector(m), p_basis_integrated_vector(n);
    auto f_functional = make_vector_functional(velocity_space, f_vector);
    f_functional.append(LocalElementIntegralFunctional<E, d>(
        local_binary_to_unary_element_integrand(LocalProductIntegrand<E, d>(), problem_.rhs_f())));
    A_operator.append(f_functional);
    auto p_basis_integrated_functional = make_vector_functional(pressure_space, p_basis_integrated_vector);
    XT::Functions::ConstantGridFunction<E> one_function(1);
    p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
        local_binary_to_unary_element_integrand(LocalProductIntegrand<E, 1>(), one_function)));
    B_operator.append(p_basis_integrated_functional);
    // Dirichlet constrainst for u
    DirichletConstraints<I, VelocitySpace> dirichlet_constraints(problem_.boundary_info(), velocity_space);
    A_operator.append(dirichlet_constraints);
    // assemble everything
    A_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    EXPECT_TRUE(is_symmetric(A));
    B_operator.assemble(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    Vector dirichlet_vector(m, 0.), reference_solution_u_vector(m, 0.), reference_solution_p_vector(n, 0.);
    auto discrete_dirichlet_values = make_discrete_function(velocity_space, dirichlet_vector);
    auto reference_solution_u = make_discrete_function(velocity_space, reference_solution_u_vector);
    auto reference_solution_p = make_discrete_function(pressure_space, reference_solution_p_vector);
    boundary_interpolation(
        problem_.dirichlet(), discrete_dirichlet_values, problem_.boundary_info(), XT::Grid::DirichletBoundary());
    default_interpolation(problem_.reference_solution_u(), reference_solution_u);
    default_interpolation(problem_.reference_solution_p(), reference_solution_p);
    Vector rhs_vector_u(m), rhs_vector_p(n);
    // create rhs vector f - A g_D for u
    A.mv(dirichlet_vector, rhs_vector_u);
    rhs_vector_u *= -1.;
    rhs_vector_u += f_vector;
    // create rhs vector -B^T g_D for p
    B.mtv(dirichlet_vector, rhs_vector_p);
    rhs_vector_p *= -1;
    // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
    // Dirichlet DoF, so we also need to clear the row of B.
    dirichlet_constraints.apply(A, rhs_vector_u);
    EXPECT_TRUE(is_symmetric(A));
    for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs())
      B.clear_row(DoF);
    EXPECT_TRUE(is_positive_definite(A));

    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
    // [A B; B^T 0] to the unit vector.
    Matrix C(n, n, pattern_C);
    size_t dof_index = 0;
    B.clear_col(dof_index);
    rhs_vector_p.set_entry(dof_index, 0.);
    C.set_entry(dof_index, dof_index, 1.);

    // now solve the system
    XT::LA::SaddlePointSolver<Vector, Matrix> solver(A, B, B, C);
    Vector solution_u(m), solution_p(n);
    // solve both by direct solver and by schurcomplement (where the schur complement is inverted by CG and the inner
    // solves with A are using a direct method)
    for (std::string type : {"direct", "cg_direct_schurcomplement"}) {
      solver.apply(rhs_vector_u, rhs_vector_p, solution_u, solution_p, type);

      // add dirichlet values to u
      const auto actual_u_vector = solution_u + dirichlet_vector;
      // ensure int_\Omega p = 0
      auto p_integral = p_basis_integrated_vector * solution_p;
      auto p_ref_integral = p_basis_integrated_vector * reference_solution_p_vector;
      auto p_correction = p_basis_integrated_vector;
      auto p_correction_func = make_discrete_function(pressure_space, p_correction);
      auto p_ref_correction = reference_solution_p_vector;
      auto p_ref_correction_func = make_discrete_function(pressure_space, p_ref_correction);
      const auto vol_domain = 4.;
      XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain);
      XT::Functions::ConstantGridFunction<E> const_p_ref_integral_func(p_ref_integral / vol_domain);
      default_interpolation(const_p_integral_func, p_correction_func);
      default_interpolation(const_p_ref_integral_func, p_ref_correction_func);
      const auto actual_p_vector = solution_p - p_correction;
      const auto actual_p_ref_vector = reference_solution_p_vector - p_ref_correction;
      // calculate difference to reference solution
      const auto u_diff_vector = actual_u_vector - reference_solution_u_vector;
      const auto p_diff_vector = actual_p_vector - actual_p_ref_vector;

      auto sol_u_func = make_discrete_function(velocity_space, actual_u_vector);
      auto sol_p_func = make_discrete_function(pressure_space, actual_p_vector);
      auto p_diff = make_discrete_function(pressure_space, p_diff_vector);
      auto u_diff = make_discrete_function(velocity_space, u_diff_vector);
      auto actual_p_ref = make_discrete_function(pressure_space, actual_p_ref_vector);
      bool visualize = true;
      std::string grid_name = XT::Common::Typename<G>::value();
      if (visualize) {
        sol_u_func.visualize("solution_u_" + type + "_" + grid_name);
        sol_p_func.visualize("solution_p_" + type + "_" + grid_name);
        reference_solution_u.visualize("u_ref_" + grid_name);
        actual_p_ref.visualize("p_ref_" + grid_name);
        u_diff.visualize("u_error_" + type + "_" + grid_name);
        p_diff.visualize("p_error_" + type + "_" + grid_name);
      }
      DXTC_EXPECT_FLOAT_LE(l2_norm(problem_.grid_view(), u_diff), expected_error_u);
      DXTC_EXPECT_FLOAT_LE(l2_norm(problem_.grid_view(), p_diff), expected_error_p);
    }
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
    : BaseType(ProblemType(ProblemType::default_diffusion(),
                           ProblemType::default_rhs_f(),
                           ProblemType::default_rhs_g(),
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
