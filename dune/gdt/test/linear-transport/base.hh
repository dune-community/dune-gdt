// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TEST_LINEAR_TRANSPORT_BASE_HH
#define DUNE_GDT_TEST_LINEAR_TRANSPORT_BASE_HH

#include <cmath>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/fvector.hh>
#include <dune/xt/functions/lambda/function.hh>
#include <dune/xt/grid/gridprovider/cube.hh>

#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Test {


template <class G>
struct LinearTransportProblem
{
  static const constexpr size_t d = G::dimension;
  using DomainType = XT::Common::FieldVector<double, d>;

  const DomainType direction;
  const XT::Functions::LambdaFunction<1, d, 1> flux;
  const double T_end;

  LinearTransportProblem()
    : direction(XT::Common::from_string<DomainType>("[1 0 0]"))
    , flux(1,
           [&](const auto& u, const auto& /*param*/) { return direction * u; },
           "linear_transport",
           {},
           [&](const auto& /*u*/, const auto& /*param*/) { return direction; })
    , T_end(1.)
  {
  }

  XT::Grid::GridProvider<G> make_initial_grid() const
  {
    return XT::Grid::make_cube_grid<G>(0., 1., 16u);
  }

  template <class Vector, class GV>
  DiscreteFunction<Vector, GV> make_exact_solution__periodic_boundaries(const SpaceInterface<GV, 1>& space,
                                                                        const double& time) const
  {
    DUNE_THROW_IF(time > 1., XT::Common::Exceptions::wrong_input_given, "time = " << time);
    const auto indicator = [](const auto& x) {
      if (0.25 <= x && x <= 0.5)
        return 1.;
      else
        return 0.;
    };
    return interpolate<Vector>(
        0, [&](const auto& xx, const auto& /*mu*/) { return indicator(std::fmod(xx[0] - time + 1., 1.)); }, space);
  } // ... make_exact_solution__periodic_boundaries(...)
}; // struct LinearTransportProblem


} // namespace Test
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_LINEAR_TRANSPORT_BASE_HH
