// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH
#define DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH

#include <cmath>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/functions/lambda/function.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/tools/euler.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct InviscidCompressibleFlowEulerProblem
{
  static const constexpr size_t d = G::dimension;
  static const constexpr size_t m = d + 2;
  using DomainType = XT::Common::FieldVector<double, d>;
  using RangeType = XT::Common::FieldVector<double, m>;

  const EulerTools<d> euler_tools;
  const XT::Functions::LambdaFunction<m, d, m> flux;
  const double T_end;

  InviscidCompressibleFlowEulerProblem()
    : euler_tools(1.4) // air or water at roughly 20 deg Cels.
    , flux(euler_tools.flux_order(),
           [&](const auto& u, const auto& /*param*/) { return euler_tools.flux(u); },
           "euler_flux",
           {},
           [&](const auto& u, const auto& /*param*/) { return euler_tools.flux_jacobian(u); })
    , T_end(1.)
  {
  }

  XT::Grid::GridProvider<G> make_initial_grid() const
  {
    return XT::Grid::make_cube_grid<G>(-1., 1., 16u);
  }

  template <class Vector, class GV>
  DiscreteFunction<Vector, GV, m> make_initial_values(const SpaceInterface<GV, m>& space) const
  {
    return interpolate<Vector>(
        0,
        [&](const auto& xx, const auto& /*mu*/) {
          RangeType primitive_variables(0.);
          // density
          if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
            primitive_variables[0] = 4.;
          else
            primitive_variables[0] = 1.;
          // velocity
          for (size_t ii = 0; ii < d; ++ii)
            primitive_variables[1 + ii] = 0.;
          // pressure
          if (XT::Common::FloatCmp::ge(xx, DomainType(-0.5)) && XT::Common::FloatCmp::le(xx, DomainType(0)))
            primitive_variables[m - 1] = 1.6;
          else
            primitive_variables[m - 1] = 0.4;
          return euler_tools.to_conservative(primitive_variables);
        },
        space);
  } // ... make_initial_values(...)
}; // struct InviscidCompressibleFlowEulerProblem


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_INVISCID_COMPRESSIBLE_FLOW_BASE_HH
