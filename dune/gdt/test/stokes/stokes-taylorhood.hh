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

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/base/function-as-grid-function.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/test/stationary-eocstudies/diffusion-ipdg.hh>

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
      XT::Functions::FunctionInterface<d, 1, 1> diffusion_factor,
      XT::Functions::FunctionInterface<d, d, 1> rhs = XT::Functions::ConstantFunction<d, d, 1>(0., "zero rhs"),
      XT::Functions::FunctionInterface<d, d, 1> dirichlet =
          XT::Functions::ConstantFunction<d, d, 1>(0., "dirichlet zero boundary values"))
    : diffusion_factor_(XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, RangeField>(diffusion_factor))
    , rhs_(XT::Functions::FunctionAsGridFunctionWrapper<E, d, 1, RangeField>(rhs))
    , dirichlet_(XT::Functions::FunctionAsGridFunctionWrapper<E, d, 1, RangeField>(dirichlet))
    , boundary_info_()
  {}

  XT::Grid::GridProvider<G> make_initial_grid()
  {
    if (std::is_same<G, YASP_2D_EQUIDISTANT_OFFSET>::value) {
      return XT::Grid::make_cube_grid<G>(-1, 1, 8);
#if HAVE_DUNE_ALUGRID
    } else if (std::is_same<G, ALU_2D_SIMPLEX_CONFORMING>::value) {
      auto grid = XT::Grid::make_cube_grid<G>(-1, 1, 4);
      grid.global_refine(2);
      return grid;
    } else if (std::is_same<G, ALU_2D_SIMPLEX_NONCONFORMING>::value) {
      return XT::Grid::make_cube_grid<G>(-1, 1, 8);
    } else if (std::is_same<G, ALU_2D_CUBE>::value) {
      auto grid = XT::Grid::make_cube_grid<G>(-1, 1, 8);
      grid.global_refine(1);
      return grid;
#endif // HAVE_DUNE_ALUGRID
    } else
      EXPECT_TRUE(false) << "Please add a specialization for '" << XT::Common::Typename<G>::value << "'!";
  } // ... make_initial_grid(...)

  const XT::Functions::FunctionAsGridFunctionWrapper<E, 1, 1, RangeField> diffusion_factor_;
  const XT::Functions::FunctionAsGridFunctionWrapper<E, d, 1, RangeField> rhs_;
  const XT::Functions::FunctionAsGridFunctionWrapper<E, d, 1, RangeField> dirichlet_;
  const XT::Grid::AllDirichletBoundaryInfo<I> boundary_info_;
}; // class ESV2007DiffusionProblem


template <class G>
class StokesDirichletTest : public StationaryEocStudy<G>
{
  using BaseType = StationaryEocStudy<G>;

  using BaseType::d;
  using typename BaseType::E;
  using typename BaseType::FF;
  using typename BaseType::FT;
  using typename BaseType::GP;
  using typename BaseType::GV;
  using typename BaseType::I;
  using typename BaseType::V;

public:
  ESV2007DiffusionTest()
    : BaseType()
    , problem()
  {}

protected:
  std::vector<std::string> norms() const override final
  {
    return {"H_1_semi"};
  }

  void compute_reference_solution() override final
  {
    auto& self = *this;
    if (self.reference_solution_on_reference_grid_)
      return;
    self.reference_grid_ = std::make_unique<GP>(make_initial_grid());
    for (size_t ref = 0; ref < self.num_refinements_ + self.num_additional_refinements_for_reference_; ++ref)
      self.reference_grid_->global_refine(DGFGridInfo<G>::refineStepsForHalf());
    auto backup_space_type = self.space_type_;
    self.space_type_ = "dg_p" + DXTC_TEST_CONFIG_GET("setup.reference_solution_order", "3");
    self.reference_space_ = self.make_space(*self.reference_grid_);
    self.space_type_ = backup_space_type;
    self.reference_solution_on_reference_grid_ = std::make_unique<V>(
        interpolate<V>(XT::Functions::ESV2007::Testcase1ExactSolution<d, 1>(), *self.reference_space_).dofs().vector());
    // visualize
    self.visualize_(
        make_discrete_function(*self.reference_space_, *self.reference_solution_on_reference_grid_),
        "reference_solution_on_refinement_"
            + XT::Common::to_string(self.num_refinements_ + self.num_additional_refinements_for_reference_));
  } // ... compute_reference_solution(...)

  const XT::Grid::BoundaryInfo<I>& boundary_info() const override final
  {
    return problem.boundary_info;
  }

  const FF& diffusion_factor() const override final
  {
    return problem.diffusion_factor.template as_grid_function<E>();
  }

  const FT& diffusion_tensor() const override final
  {
    return problem.diffusion_tensor.template as_grid_function<E>();
  }

  const FF& force() const override final
  {
    return problem.force.template as_grid_function<E>();
  }

  GP make_initial_grid() override final
  {
    return problem.make_initial_grid();
  }

  ESV2007DiffusionProblem<GV> problem;
}; // struct ESV2007DiffusionTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_HEAT_EQUATION_ESV2007_HH
