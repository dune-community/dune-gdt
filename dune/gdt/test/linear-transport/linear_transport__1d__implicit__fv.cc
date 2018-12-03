// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_TIMED_LOGGING 1
#define DUNE_XT_COMMON_TEST_MAIN_ENABLE_INFO_LOGGING 1

#include <dune/xt/common/test/main.hxx> // <- this one has to come first (includes the config.h)!

#include <dune/xt/grid/grids.hh>

#include <dune/gdt/test/instationary-eocstudies/hyperbolic-nonconforming.hh>

#include "base.hh"

using namespace Dune;
using namespace Dune::GDT::Test;


// This is just for comparison
using LinearTransport1dExplicitFvTest = LinearTransportExplicitTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(LinearTransport1dExplicitFvTest, dt_equals_h)
{
  this->visualization_steps_ =
      DXTC_CONFIG_GET("visualization_steps__explicit__dt_equals_h", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  this->run();
}


using LinearTransport1dImplicitFvTest = LinearTransportImplicitTest<YASP_1D_EQUIDISTANT_OFFSET>;
TEST_F(LinearTransport1dImplicitFvTest, dt_equals_h)
{
  this->visualization_steps_ =
      DXTC_CONFIG_GET("visualization_steps__implicit__dt_equals_h", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  this->run();
}


// This one is also just for comparison
class LinearTransport1dExplicitWithAutomaticDtTest : public LinearTransportExplicitTest<YASP_1D_EQUIDISTANT_OFFSET>
{
protected:
  using BaseBaseType = InstationaryNonconformingHyperbolicEocStudy<YASP_1D_EQUIDISTANT_OFFSET, 1>;
  using typename BaseBaseType::S;

  XT::LA::ListVectorArray<V> solve(const S& space, const double T_end) override final
  {
    const auto u_0 = this->make_initial_values(space);
    const auto op = this->make_lhs_operator(space);
    const auto dt = this->estimate_fixed_explicit_fv_dt(space);
    this->current_data_["quantity"]["dt"] = dt;
    this->current_data_["quantity"]["explicit_fv_dt"] = dt;
    return solve_instationary_system_explicit_euler(u_0, *op, T_end, dt);
  }
}; // class LinearTransport1dExplicitWithAutomaticDtTest

TEST_F(LinearTransport1dExplicitWithAutomaticDtTest, automatic_dt)
{
  this->visualization_steps_ =
      DXTC_CONFIG_GET("visualization_steps__explicit__automatic_dt", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  this->run();
}


TEST_F(LinearTransport1dImplicitFvTest, larger_dt)
{
  this->dt_factor_ = DXTC_CONFIG_GET("implicit_dt_factor", 10.);
  this->visualization_steps_ =
      DXTC_CONFIG_GET("visualization_steps__implicit__larger_dt", 0); // <- Something like 100 to visualize!
  this->space_type_ = "fv";
  this->numerical_flux_type_ = "upwind";
  this->run();
}
