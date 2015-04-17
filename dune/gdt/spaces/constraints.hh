// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_HH

#include <mutex>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/la/container/interfaces.hh>

namespace Dune {
namespace GDT {
namespace internal {


// forward, needed for friendlyness
template< class TestSpaceType, class AnsatzSpaceType, class GridViewType, class ConstraintsType >
class ConstraintsWrapper;


} // namespace internal
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
  typedef DirichletConstraints< IntersectionType >                 ThisType;
public:
  typedef internal::DirichletConstraintsTraits< IntersectionType > Traits;
  typedef Stuff::Grid::BoundaryInfoInterface< IntersectionType >   BoundaryInfoType;

  DirichletConstraints(const BoundaryInfoType& bnd_info, const size_t sz, const bool set = true)
    : boundary_info_(bnd_info)
    , size_(sz)
    , set_(set)
  {}

  // manual copy ctor needed bc. of the mutex
  DirichletConstraints(const ThisType& other)
    : boundary_info_(other.boundary_info_)
    , size_(other.size_)
    , set_(other.set_)
    , dirichlet_DoFs_(other.dirichlet_DoFs_)
  {}

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  size_t size() const
  {
    return size_;
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
  template< class T, class A, class GV, class C >
  friend class GDT::internal::ConstraintsWrapper;

  const BoundaryInfoType& boundary_info_;
  const size_t size_;
  const bool set_;
  std::set< size_t > dirichlet_DoFs_;
  std::mutex mutex_;
}; // class DirichletConstraints


} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
