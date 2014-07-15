// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_ASSEMBLER_FUNCTORS_HH
#define DUNE_GDT_ASSEMBLER_FUNCTORS_HH

#include <dune/stuff/grid/boundaryinfo.hh>

#include "../../assembler/gridwalker.hh"

namespace Dune {
namespace GDT {
namespace Functor {


template< class GridViewImp >
class DirichletDetector
  : public Codim1< GridViewImp >
{
  typedef Codim1< GridViewImp > BaseType;
public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;
  typedef typename BaseType::IntersectionType IntersectionType;

  DirichletDetector(const Stuff::Grid::BoundaryInfoInterface< IntersectionType >& boundary_info)
    : boundary_info_(boundary_info)
    , found_(0)
  {}

  virtual ~DirichletDetector() {}

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& /*inside_entity*/,
                           const EntityType& /*outside_entity*/) DS_OVERRIDE
  {
    if (boundary_info_.dirichlet(intersection))
      ++found_;
  }

  bool found() const
  {
    return found_ > 0;
  }

private:
  const Stuff::Grid::BoundaryInfoInterface< IntersectionType >& boundary_info_;
  size_t found_;
}; // class DirichletDetector


} // namespace Functor
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_FUNCTORS_HH
