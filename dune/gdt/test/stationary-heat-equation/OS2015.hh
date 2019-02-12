// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_STATIONARY_HEAT_EQUATION_OS2015_HH
#define DUNE_GDT_TEST_STATIONARY_HEAT_EQUATION_OS2015_HH

#include <dune/xt/common/test/gtest/gtest.h>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/eye-matrix.hh>
#include <dune/xt/grid/boundaryinfo/alldirichlet.hh>
#include <dune/xt/grid/gridprovider/cube.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/grids.hh>
#include <dune/xt/functions/constant.hh>
#include <dune/xt/functions/indicator.hh>
#include <dune/xt/functions/spe10/model1.hh>

#include <dune/gdt/interpolations/default.hh>
#include <dune/gdt/test/stationary-eocstudies/diffusion-ipdg.hh>

namespace Dune {
namespace GDT {
namespace Test {


/**
 * \brief Problem definition from [1], Section 6.1
 *
 * [1]: Ohlberger, Schindler, 2015, Error Control for the Localized Reduced Basis Multiscale Method with Adaptive
 *      On-Line Enrichment, SIAM J. Sci. Comput., Vol. 37, No. 6, pp. A2865-A2895
 */
template <class GV>
struct OS2015MultiscaleProblem
{
  static_assert(XT::Grid::is_view<GV>::value, "");
  static_assert(GV::dimension == 2, "");

  static const constexpr size_t d = GV::dimension;
  using E = XT::Grid::extract_entity_t<GV>;
  using I = XT::Grid::extract_intersection_t<GV>;
  using G = typename GV::Grid;

  OS2015MultiscaleProblem()
    : diffusion_factor(1)
    , diffusion_tensor(XT::Data::spe10_model1_filename(), {0., 0.}, {5., 1.})
    , dirichlet(0)
    , neumann(0)
    , force({/*domain_with_positive_value=*/{{{0.95, 1.1}, {0.3, 0.45}}, 2e3},
             /*first_domain_with_negative_value=*/{{{3., 3.15}, {0.75, 0.9}}, -1e3},
             /*second_domain_with_negative_value=*/{{{4.25, 4.4}, {0.25, 0.4}}, -1e3}})
  {}

  XT::Grid::GridProvider<G> make_initial_grid()
  {
    if (std::is_same<G, YASP_2D_EQUIDISTANT_OFFSET>::value) {
      return XT::Grid::make_cube_grid<G>({0., 0.}, {5., 1.}, {100u, 20u});
#if HAVE_DUNE_ALUGRID
    } else if (std::is_same<G, ALU_2D_SIMPLEX_CONFORMING>::value) {
      auto grid = XT::Grid::make_cube_grid<G>({0., 0.}, {5., 1.}, {100u, 20u});
      grid.global_refine(2);
      return grid;
#endif // HAVE_DUNE_ALUGRID
    } else
      EXPECT_TRUE(false) << "Please add a specialization for '" << XT::Common::Typename<G>::value << "'!";
  } // ... make_initial_grid(...)

  const XT::Functions::ConstantFunction<d> diffusion_factor;
  const XT::Functions::Spe10::Model1Function<E, d, d> diffusion_tensor;
  const XT::Functions::ConstantFunction<d> dirichlet;
  const XT::Functions::ConstantFunction<d> neumann;
  const XT::Functions::IndicatorFunction<d> force;
  const XT::Grid::AllDirichletBoundaryInfo<I> boundary_info;
}; // class OS2015MultiscaleProblem


/**
 * \brief To reproduce the results in [1], Section 6.1
 *
 * [1]: Ohlberger, Schindler, 2015, Error Control for the Localized Reduced Basis Multiscale Method with Adaptive
 *      On-Line Enrichment, SIAM J. Sci. Comput., Vol. 37, No. 6, pp. A2865-A2895
 */
template <class G>
class OS2015MultiscaleTest : public StationaryDiffusionIpdgEocStudy<G>
{
  using BaseType = StationaryDiffusionIpdgEocStudy<G>;

  using BaseType::d;
  using typename BaseType::E;
  using typename BaseType::FF;
  using typename BaseType::FT;
  using typename BaseType::GP;
  using typename BaseType::GV;
  using typename BaseType::I;
  using typename BaseType::V;

public:
  OS2015MultiscaleTest()
    : BaseType()
    , problem()
  {}

protected:
  std::vector<std::string> norms() const override final
  {
    return {"eta_NC", "eta_R", "eta_DF"};
  }

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
    return problem.diffusion_tensor;
  }

  const FF& force() const override final
  {
    return problem.force.template as_grid_function<E>();
  }

  GP make_initial_grid() override final
  {
    return problem.make_initial_grid();
  }

  OS2015MultiscaleProblem<GV> problem;
}; // struct OS2015MultiscaleTest


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_STATIONARY_HEAT_EQUATION_OS2015_HH
