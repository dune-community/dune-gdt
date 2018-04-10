// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2017)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_PROLONGATIONS_HH
#define DUNE_GDT_PROLONGATIONS_HH

#include <dune/xt/functions/reinterpret.hh>
#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/discretefunction/reinterpret.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [most general
 *        variant].
 *
 * \note This does not clear target.dofs().vector(). Thus, if prolongation_grid_view only covers a part of the domain of
 *       target.space().grid_view(), other contributions in target remain (which is on purpose).
 */
template <class SV, class SGV, size_t r, size_t rC, class SR, class TV, class TGV, class TR, class PGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<TGV>, typename PGV::Grid::template Codim<0>::Entity>::value,
                 void>
prolong(const DiscreteFunction<SV, SGV, r, rC, SR>& source,
        DiscreteFunction<TV, TGV, r, rC, TR>& target,
        const GridView<PGV>& prolongation_grid_view)
{
  interpolate(reinterpret(source, prolongation_grid_view), target, prolongation_grid_view);
}


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [uses
 *        target.space().grid_view() as prolongation_grid_view].
 */
template <class SV, class SGV, size_t r, size_t rC, class SR, class TV, class TGV, class TR>
void prolong(const DiscreteFunction<SV, SGV, r, rC, SR>& source, DiscreteFunction<TV, TGV, r, rC, TR>& target)
{
  prolong(source, target, target.space().grid_view());
}


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [creates
 *        suitable target_function, TargetVectorType has to be provided].
 *
 * Use as in
\code
auto target_function = prolong<TargetVectorType>(source, target_space, prolongation_grid_view);
\endcode
 */
template <class TargetVectorType, class SV, class SGV, size_t r, size_t rC, class SR, class TGV, class TR, class PGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<TGV>, typename PGV::Grid::template Codim<0>::Entity>::value,
                 DiscreteFunction<TargetVectorType, TGV, r, rC, TR>>
prolong(const DiscreteFunction<SV, SGV, r, rC, SR>& source,
        const SpaceInterface<TGV, r, rC, TR>& target_space,
        const GridView<PGV>& prolongation_grid_view)
{
  auto target_function = make_discrete_function<TargetVectorType>(target_space);
  prolong(source, target_function, prolongation_grid_view);
  return target_function;
}


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [creates
 *        suitable target_function, TargetVectorType has to be provided, uses target.space().grid_view() as
 *        prolongation_grid_view].
 *
 * Use as in
\code
auto target_function = prolong<TargetVectorType>(source, target_space);
\endcode
 */
template <class TargetVectorType, class SV, class SGV, size_t r, size_t rC, class SR, class TGV, class TR>
DiscreteFunction<TargetVectorType, TGV, r, rC, TR> prolong(const DiscreteFunction<SV, SGV, r, rC, SR>& source,
                                                           const SpaceInterface<TGV, r, rC, TR>& target_space)
{
  auto target_function = make_discrete_function<TargetVectorType>(target_space);
  prolong(source, target_function);
  return target_function;
}


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [creates
 *        suitable target_function with same VectorType as source].
 */
template <class V, class SGV, size_t r, size_t rC, class SR, class TGV, class TR, class PGV>
std::enable_if_t<std::is_same<XT::Grid::extract_entity_t<TGV>, typename PGV::Grid::template Codim<0>::Entity>::value,
                 DiscreteFunction<V, TGV, r, rC, TR>>
prolong(const DiscreteFunction<V, SGV, r, rC, SR>& source,
        const SpaceInterface<TGV, r, rC, TR>& target_space,
        const GridView<PGV>& prolongation_grid_view)
{
  auto target_function = make_discrete_function<V>(target_space);
  prolong(source, target_function, prolongation_grid_view);
  return target_function;
}


/**
 * \brief Prolongs a DiscreteFunction from one (usually coarser) GridView onto another (usually finer) one [creates
 *        suitable target_function with same VectorType as source, uses target.space().grid_view() as
 *        prolongation_grid_view].
 */
template <class V, class SGV, size_t r, size_t rC, class SR, class TGV, class TR>
DiscreteFunction<V, TGV, r, rC, TR> prolong(const DiscreteFunction<V, SGV, r, rC, SR>& source,
                                            const SpaceInterface<TGV, r, rC, TR>& target_space)
{
  auto target_function = make_discrete_function<V>(target_space);
  prolong(source, target_function);
  return target_function;
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PROLONGATIONS_HH
