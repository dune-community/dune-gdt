// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   René Fritze     (2016, 2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
#define DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/base/reinterpret.hh>

#include "default.hh"

namespace Dune {
namespace GDT {


/**
 * \sa XT::Functions::ReinterpretLocalizableFunction
 * \sa XT::Functions::reinterpret
 */
template <class TargetElement, class SGV, size_t r, size_t rC, class R, class V>
XT::Functions::ReinterpretLocalizableFunction<SGV, TargetElement, r, rC, R>
reinterpret(const DiscreteFunction<V, SGV, r, rC, R>& source)
{
  return XT::Functions::ReinterpretLocalizableFunction<SGV, TargetElement, r, rC, R>(source,
                                                                                     source.space().grid_view());
}


/**
 * \sa XT::Functions::ReinterpretLocalizableFunction
 * \sa XT::Functions::reinterpret
 */
template <class SGV, size_t r, size_t rC, class R, class V, class TargetGridView>
std::enable_if_t<
    XT::Grid::is_layer<TargetGridView>::value,
    XT::Functions::ReinterpretLocalizableFunction<SGV, XT::Grid::extract_entity_t<TargetGridView>, r, rC, R>>
reinterpret(const DiscreteFunction<V, SGV, r, rC, R>& source, const TargetGridView& /*target_grid_view*/)
{
  return XT::Functions::ReinterpretLocalizableFunction<SGV, XT::Grid::extract_entity_t<TargetGridView>, r, rC, R>(
      source, source.space().grid_view());
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
