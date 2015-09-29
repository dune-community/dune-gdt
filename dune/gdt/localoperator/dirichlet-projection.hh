// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALOPERATOR_DIRICHLET_PROJECTION_HH
#define DUNE_GDT_LOCALOPERATOR_DIRICHLET_PROJECTION_HH

#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/discretefunction/local.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward
template< class BoundaryInfoType >
class LocalDirichletProjectionOperator;


namespace internal {


/**
 * \todo Check BoundaryInfoType!
 */
template< class BoundaryInfoType >
class LocalDirichletProjectionOperatorTraits
{
public:
  typedef LocalDirichletProjectionOperator< BoundaryInfoType > derived_type;
};


} // namespace internal


template< class BoundaryInfoType >
class LocalDirichletProjectionOperator
    : public LocalOperatorInterface< internal::LocalDirichletProjectionOperatorTraits< BoundaryInfoType > >
{
public:
  typedef internal::LocalDirichletProjectionOperatorTraits< BoundaryInfoType > Traits;

  LocalDirichletProjectionOperator(const BoundaryInfoType& boundary_info)
    : boundary_info_(boundary_info)
  {}

  template< class E, class D, size_t d, class R, size_t r, size_t rC, class RS, class V >
  void apply(const Stuff::LocalizableFunctionInterface< E, D, d, R, r, rC >& source,
             LocalDiscreteFunction< RS, V >& local_range) const
  {
    const auto& entity = local_range.entity();
    if (!entity.hasBoundaryIntersections())
      return;
    const auto local_dirichlet_DoFs = local_range.space().local_dirichlet_DoFs(entity, boundary_info_);
    if (local_dirichlet_DoFs.size() == 0)
      return;
    const auto lagrange_points = local_range.space().lagrange_points(entity);
    const auto local_source = source.local_function(entity);
    auto& local_range_DoF_vector = local_range.vector();
    assert(lagrange_points.size() == local_range_DoF_vector.size());
    for (const size_t& local_DoF_id : local_dirichlet_DoFs)
      local_range_DoF_vector.set(local_DoF_id, local_source->evaluate(lagrange_points[local_DoF_id]));
  } // ... apply(...)

private:
  const BoundaryInfoType& boundary_info_;
}; // class LocalDirichletProjectionOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALOPERATOR_DIRICHLET_PROJECTION_HH
