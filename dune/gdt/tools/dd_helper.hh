// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_TOOLS_DD_HELPER
#define DUNE_GDT_TOOLS_DD_HELPER

#include <memory>

#include <dune/common/dynvector.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container/conversion.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/h1/continuous-lagrange.hh>
#include <dune/gdt/spaces/l2/discontinuous-lagrange.hh>

namespace Dune {
namespace GDT {

template <class GV>
std::unique_ptr<GDT::SpaceInterface<GV>> make_subdomain_space(GV subdomain_grid_view, const std::string& space_type)
{
  if (space_type.size() >= 4 && space_type.substr(0, 4) == "cg_p") {
    const auto order = XT::Common::from_string<int>(space_type.substr(4));
    return std::make_unique<ContinuousLagrangeSpace<GV>>(subdomain_grid_view, order);
  } else if (space_type.size() >= 4 && space_type.substr(0, 4) == "dg_p") {
    const auto order = XT::Common::from_string<int>(space_type.substr(4));
    return std::make_unique<DiscontinuousLagrangeSpace<GV>>(subdomain_grid_view, order);
  } else
    DUNE_THROW(XT::Common::Exceptions::wrong_input_given,
               "space_type = " << space_type << "\n   has to be 'cg_pX' or 'dg_pX' for some order X!");
  return nullptr;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_DD_HELPER
