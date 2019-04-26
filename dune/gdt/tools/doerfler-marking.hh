// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_TOOLS_DOERFLER_MARKING_HH
#define DUNE_GDT_TOOLS_DOERFLER_MARKING_HH

#include <algorithm>
#include <set>
#include <utility>
#include <vector>

#include <dune/xt/la/container/vector-interface.hh>

namespace Dune {
namespace GDT {


template <class V>
std::pair<std::set<size_t>, std::set<size_t>>
doerfler_marking(const XT::LA::VectorInterface<V>& squared_local_indicators,
                 const double& refine_factor,
                 const double& coarsen_factor = 0.)
{
  // sort
  std::vector<std::pair<typename XT::LA::VectorInterface<V>::ScalarType, size_t>> accumulated_sorted_local_indicators(
      squared_local_indicators.size());
  for (size_t ii = 0; ii < squared_local_indicators.size(); ++ii)
    accumulated_sorted_local_indicators[ii] = {std::pow(squared_local_indicators[ii], 2), ii};
  std::sort(accumulated_sorted_local_indicators.begin(),
            accumulated_sorted_local_indicators.end(),
            [](const auto& a, const auto& b) { return a.first < b.first; });
  for (size_t ii = 1; ii < accumulated_sorted_local_indicators.size(); ++ii)
    accumulated_sorted_local_indicators[ii].first += accumulated_sorted_local_indicators[ii - 1].first;
  // select smallest coarsen_factor contributions for coarsening
  std::set<size_t> elements_to_be_coarsened;
  const double total_indicators = std::pow(squared_local_indicators.l2_norm(), 2);
  for (const auto& indicator : accumulated_sorted_local_indicators)
    if (indicator.first < (coarsen_factor * total_indicators))
      elements_to_be_coarsened.insert(indicator.second);
    else
      break;
  // select largest refine_factor contributions
  std::set<size_t> elements_to_be_refined;
  for (const auto& indicator : accumulated_sorted_local_indicators)
    if (indicator.first > ((1 - refine_factor) * total_indicators))
      elements_to_be_refined.insert(indicator.second);
  return {elements_to_be_coarsened, elements_to_be_refined};
} // ... doerfler_marking(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TOOLS_DOERFLER_MARKING_HH
