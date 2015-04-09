// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_HH

#include <ostream>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/interfaces.hh>

namespace Dune {
namespace GDT {
namespace Spaces {


/**
 * \brief CRTP interface for all implementations of constraints.
 *
 *        We need this interface for template matching in the SystemAssembler.
 */
template< class Traits >
class ConstraintsInterface
  : public Stuff::CRTPInterface< ConstraintsInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
}; // class ConstraintsInterface


// forward
template< class IntersectionType >
class DirichletConstraints;


namespace internal {


template< class IntersectionType >
class DirichletConstraintsTraits
{
public:
  typedef DirichletConstraints< IntersectionType > derived_type;
};


} // namespace internal


template< class IntersectionType >
class DirichletConstraints
  : public ConstraintsInterface< internal::DirichletConstraintsTraits< IntersectionType > >
{
public:
  typedef internal::DirichletConstraintsTraits< IntersectionType > Traits;
  typedef Stuff::Grid::BoundaryInfoInterface< IntersectionType >   BoundaryInfoType;

  DirichletConstraints(const BoundaryInfoType& bnd_info, const size_t sz, const bool set = true)
    : boundary_info_(bnd_info)
    , size_(sz)
    , set_(set)
  {}

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  inline void insert(const size_t DoF)
  {
    assert(DoF < size_);
    dirichlet_DoFs_.insert(DoF);
  }

  template< class M >
  void apply(Stuff::LA::MatrixInterface< M >& matrix) const
  {
    assert(matrix.rows() == size_);
    if (set_) {
      for (const auto& DoF : dirichlet_DoFs_)
        matrix.unit_row(DoF);
    } else {
      for (const auto& DoF : dirichlet_DoFs_)
        matrix.clear_row(DoF);
    }
  } // ... apply(...)

  template< class V >
  void apply(Stuff::LA::VectorInterface< V >& vector) const
  {
    assert(vector.size() == size_);
    for (const auto& DoF : dirichlet_DoFs_)
      vector[DoF] = 0.0;
  }

  template< class M, class V >
  void apply(Stuff::LA::MatrixInterface< M >& matrix, Stuff::LA::VectorInterface< V >& vector) const
  {
    assert(matrix.rows() == size_);
    assert(vector.size() == size_);
    if (set_) {
      for (const auto& DoF : dirichlet_DoFs_) {
        matrix.unit_row(DoF);
        vector[DoF] = 0.0;
      }
    } else {
      for (const auto& DoF : dirichlet_DoFs_) {
        matrix.clear_row(DoF);
        vector[DoF] = 0.0;
      }
    }
  } // ... apply(...)

private:
  const BoundaryInfoType& boundary_info_;
  const size_t size_;
  const bool set_;
  std::set< size_t > dirichlet_DoFs_;
}; // class DirichletConstraints


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
