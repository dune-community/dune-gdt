// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2019)

#ifndef DUNE_GDT_TEST_BURGERS_BASE_HH
#define DUNE_GDT_TEST_BURGERS_BASE_HH

#include <cmath>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/common/memory.hh>
#include <dune/xt/functions/generic/function.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/test/instationary-eocstudies/hyperbolic-nonconforming.hh>

namespace Dune {
namespace GDT {


template <class G>
struct BurgersProblem
{
  static const constexpr size_t d = G::dimension;
  static_assert(d == 1, "Not implemented yet!");
  //  using DomainType = XT::Common::FieldVector<double, d>;

  const XT::Functions::GenericFunction<1, d, 1> flux;
  const double T_end;

  BurgersProblem()
    : flux(2,
           [&](const auto& u, const auto& /*param*/) { return 0.5 * u * u; },
           "burgers",
           {},
           [&](const auto& u, const auto& /*param*/) { return u; })
    , T_end(1.)
  {}

  XT::Grid::GridProvider<G> make_initial_grid() const
  {
    return XT::Grid::make_cube_grid<G>(0., 1., 16u);
  }

  template <class Vector, class GV>
  DiscreteFunction<Vector, GV> make_initial_values(const SpaceInterface<GV, 1>& space) const
  {
    std::cerr << "TODO: Use generic interpolate (?) once implemented." << std::endl;
    return default_interpolation<Vector>(3,
                                         [&](const auto& xx, const auto& /*mu*/) {
                                           return std::exp(-std::pow(xx[0] - 0.33, 2) / (2 * std::pow(0.075, 2)));
                                         },
                                         space);
  }
}; // struct BurgersProblem


template <class G>
class BurgersTest
  : XT::Common::ConstStorageProvider<BurgersProblem<G>>
  , public GDT::Test::InstationaryNonconformingHyperbolicEocStudy<G, 1>
{
  using Problem = XT::Common::ConstStorageProvider<BurgersProblem<G>>;
  using BaseType = GDT::Test::InstationaryNonconformingHyperbolicEocStudy<G, 1>;

protected:
  using BaseType::d;
  using typename BaseType::BS;
  using typename BaseType::DF;
  using typename BaseType::F;
  using typename BaseType::GP;
  using typename BaseType::S;
  using typename BaseType::V;

public:
  BurgersTest(const std::string timestepping)
    : Problem(new BurgersProblem<G>())
    , BaseType(Problem::access().T_end,
               timestepping,
               [&](const auto& solution, const auto& prefix) {
                 const auto end_time =
                     std::min(this->T_end_, this->time_points_from_vector_array(solution.dof_vectors()).back());
                 for (size_t ii = 0; ii < this->visualization_steps_; ++ii) {
                   const double time = ii * (end_time / this->visualization_steps_);
                   solution.evaluate(time).visualize(XT::Common::Test::get_unique_test_name() + "__" + prefix
                                                     + "_solution_" + XT::Common::to_string(ii));
                 }
               })
    , visualization_steps_(0)
  {}

protected:
  const F& flux() const override
  {
    return Problem::access().flux;
  }

  DF make_initial_values(const S& space) override
  {
    return Problem::access().template make_initial_values<V>(space);
  }

  GP make_initial_grid() override
  {
    return Problem::access().make_initial_grid();
  }

protected:
  size_t visualization_steps_;
}; // class BurgersTest


template <class G>
class BurgersExplicitTest : public BurgersTest<G>
{
protected:
  using BaseType = BurgersTest<G>;
  using typename BaseType::BS;
  using typename BaseType::GP;
  using typename BaseType::S;
  using typename BaseType::V;

  BurgersExplicitTest()
    : BaseType("explicit/fixed")
  {}

  virtual void compute_reference_solution() override
  {
    auto& self = *this;
    auto st = self.space_type_;
    self.space_type_ = "fv";
    BaseType::compute_reference_solution();
    self.space_type_ = st;
  }

  XT::LA::ListVectorArray<V> solve(const S& space, const double T_end) override
  {
    auto& self = *this;
    const auto fv_dt = self.estimate_fixed_explicit_fv_dt(space);
    const auto dt = self.estimate_fixed_explicit_dt_to_T_end(
        space,
        DXTC_TEST_CONFIG_GET("setup.estimate_fixed_explicit_dt.min_dt", 1e-4),
        T_end,
        DXTC_TEST_CONFIG_GET("setup.estimate_fixed_explicit_dt.max_overshoot", 1.25));
    self.current_data_["quantity"]["dt"] = dt;
    self.current_data_["quantity"]["explicit_fv_dt"] = fv_dt;
    Timer timer;
    const auto u_0 = self.make_initial_values(space);
    const auto op = self.make_lhs_operator(space);
    /// TODO: investigate why the tested dt from above does not work!
    auto solution = GDT::Test::solve_instationary_system_explicit_euler(
        u_0, *op, T_end, DXTC_TEST_CONFIG_GET("setup.dt_factor", 0.99) * dt);
    self.current_data_["quantity"]["time to solution (s)"] = timer.elapsed();
    return solution;
  } // ... solve(...)
}; // class BurgersExplicitTest


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_BURGERS_BASE_HH
