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
#include <dune/xt/grid/walker.hh>

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
#include <dune/gdt/operators/matrix.hh>
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

  static constexpr size_t d = GV::dimension;
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
  }

  const XT::Functions::GridFunctionInterface<E, d, d>& diffusion()
  {
    return *diffusion_;
  }

  const VectorGridFunction& rhs_f()
  {
    return *rhs_f_;
  }

  const VectorGridFunction& rhs_g()
  {
    return *rhs_g_;
  }

  const VectorGridFunction& dirichlet()
  {
    return *dirichlet_;
  }

  const VectorGridFunction& reference_solution_u()
  {
    DUNE_THROW_IF(!reference_sol_u_, Dune::InvalidStateException, "No reference solution provided!");
    return *reference_sol_u_;
  }

  const ScalarGridFunction& reference_solution_p()
  {
    DUNE_THROW_IF(!reference_sol_p_, Dune::InvalidStateException, "No reference solution provided!");
    return *reference_sol_p_;
  }

  const XT::Grid::BoundaryInfo<I>& boundary_info()
  {
    return boundary_info_;
  }

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
  static constexpr size_t d = GV::dimension;

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
    // With n = velocity_space.mapper.size(), m = pressure_space.mapper().size()
    const VelocitySpace velocity_space(grid_view, velocity_order);
    const PressureSpace pressure_space(grid_view, velocity_order - 1);
    // Define bilinear forms associated with A and B
    // - calculate A_{ij} as \int \nabla v_i \nabla v_j, using associated bilinear form a
    BilinearForm<GV, d> a(grid_view);
    a += LocalElementIntegralBilinearForm<E, d>(LocalLaplaceIntegrand<E, d>(problem_.diffusion()));
    // - calculate B_{ij} as \int -\nabla p_i div(v_j), using associated bilinear form b
    BilinearForm<GV, 1, 1, d, 1> b(grid_view);
    b += LocalElementIntegralBilinearForm<E, d, 1, RangeField, RangeField, 1>(
        LocalElementAnsatzValueTestDivProductIntegrand<E>(-1.));
    // Turn these into (discrete) matrix operators (by providing discrete function spaces) which can be assembled
    auto A = a.template with<Matrix>(velocity_space, velocity_space);
    auto B = b.template with<Matrix>(pressure_space, velocity_space);
    // calculate rhs f as \int ff v
    auto f_functional = make_vector_functional<Vector>(velocity_space);
    f_functional.append(
        LocalElementIntegralFunctional<E, d>(LocalProductIntegrand<E, d>().with_ansatz(problem_.rhs_f())));
    // calculate integrated pressure space basis \int q_i
    auto p_basis_integrated_functional = make_vector_functional<Vector>(pressure_space);
    p_basis_integrated_functional.append(LocalElementIntegralFunctional<E, 1>(
        LocalProductIntegrand<E, 1>().with_ansatz(XT::Functions::ConstantGridFunction<E>(1))));
    // Dirichlet constrainst for u
    DirichletConstraints<I, VelocitySpace> dirichlet_constraints(problem_.boundary_info(), velocity_space);
    // assemble everything
    auto walker = XT::Grid::make_walker(grid_view);
    walker.append(A);
    walker.append(B);
    walker.append(f_functional);
    walker.append(p_basis_integrated_functional);
    walker.append(dirichlet_constraints);
    walker.walk(DXTC_TEST_CONFIG_GET("setup.use_tbb", true));
    EXPECT_TRUE(is_symmetric(A.matrix()));
    const auto discrete_dirichlet_values = boundary_interpolation<Vector>(
        problem_.dirichlet(), velocity_space, problem_.boundary_info(), XT::Grid::DirichletBoundary());
    auto reference_solution_u = default_interpolation<Vector>(problem_.reference_solution_u(), velocity_space);
    auto reference_solution_p = default_interpolation<Vector>(problem_.reference_solution_p(), pressure_space);
    // create rhs vector f - A g_D for u
    auto rhs_vector_u = A.apply(discrete_dirichlet_values.dofs().vector());
    rhs_vector_u *= -1.;
    rhs_vector_u += f_functional.vector();
    // create rhs vector -B^T g_D for p
    auto rhs_vector_p = B.matrix().mtv(discrete_dirichlet_values.dofs().vector()); // TODO: use B.apply_adjoint ...
    rhs_vector_p *= -1; // ... see https://zivgitlab.uni-muenster.de/ag-ohlberger/dune-community/dune-gdt/-/issues/32
    // apply dirichlet constraints for u. We need to set the whole row of (A B; B^T 0) to the unit row for each
    // Dirichlet DoF, so we also need to clear the row of B.
    dirichlet_constraints.apply(A.matrix(), rhs_vector_u);
    EXPECT_TRUE(is_symmetric(A.matrix()));
    for (const auto& DoF : dirichlet_constraints.dirichlet_DoFs())
      B.matrix().clear_row(DoF);
#if defined(NDEBUG) || HAVE_MKL || HAVE_LAPACKE
    // This check is very slow if compiled in debug mode without LAPACKE,
    // so we disable it in that case to avoid a test timeout.
    EXPECT_TRUE(is_positive_definite(A.matrix()));
#endif

    // Fix value of p at first DoF to 0 to ensure the uniqueness of the solution, i.e, we have set the m-th row of
    // [A B; B^T 0] to the unit vector.
    auto C = make_matrix_operator<Matrix>(pressure_space);
    const size_t dof_index = 0;
    B.matrix().clear_col(dof_index);
    rhs_vector_p.set_entry(dof_index, 0.);
    C.matrix().set_entry(dof_index, dof_index, 1.);

    // now solve the system
    XT::LA::SaddlePointSolver<Vector, Matrix> solver(A.matrix(), B.matrix(), B.matrix(), C.matrix());
    auto solution_u = make_discrete_function<Vector>(velocity_space);
    auto solution_p = make_discrete_function<Vector>(pressure_space);
    // solve both by direct solver and by schurcomplement (where the schur complement is inverted by CG and the inner
    // solves with A are using a direct method)
    for (std::string type : {"direct", "cg_direct_schurcomplement"}) {
      solver.apply(rhs_vector_u, rhs_vector_p, solution_u.dofs().vector(), solution_p.dofs().vector(), type);

      // add dirichlet values to u
      solution_u.dofs().vector() += discrete_dirichlet_values.dofs().vector();
      // ensure int_\Omega p = 0
      const auto p_integral = p_basis_integrated_functional.apply(solution_p);
      const auto p_ref_integral = p_basis_integrated_functional.apply(reference_solution_p);
      const auto vol_domain = 4.;
      XT::Functions::ConstantGridFunction<E> const_p_integral_func(p_integral / vol_domain);
      XT::Functions::ConstantGridFunction<E> const_p_ref_integral_func(p_ref_integral / vol_domain);
      const auto p_correction_func = default_interpolation<Vector>(const_p_integral_func, pressure_space);
      const auto p_ref_correction_func = default_interpolation<Vector>(const_p_ref_integral_func, pressure_space);
      solution_p.dofs().vector() -= p_correction_func.dofs().vector();
      reference_solution_p.dofs().vector() -= p_ref_correction_func.dofs().vector();
      // calculate difference to reference solution
      const auto p_diff = solution_p - reference_solution_p;
      const auto u_diff = solution_u - reference_solution_u;
      std::string grid_name = XT::Common::Typename<G>::value();
      if (DXTC_TEST_CONFIG_GET("visualize", true)) {
        visualize(solution_u, "solution_u_" + type + "_" + grid_name);
        visualize(solution_p, "solution_p_" + type + "_" + grid_name);
        visualize(reference_solution_u, "u_ref_" + grid_name);
        visualize(reference_solution_p, "p_ref_" + grid_name);
        visualize(u_diff, grid_view, "u_error_" + type + "_" + grid_name);
        visualize(p_diff, grid_view, "p_error_" + type + "_" + grid_name);
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
  static constexpr size_t d = ProblemType::d;

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
