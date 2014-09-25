// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_PARALLEL_HH
#define DUNE_GDT_SPACES_PARALLEL_HH

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#if HAVE_DUNE_ISTL
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif

#if HAVE_DUNE_PDELAB
#include <dune/stuff/common/disable_warnings.hh>
#include <dune/pdelab/backend/istl/parallelhelper.hh>
#include <dune/stuff/common/reenable_warnings.hh>
#endif

#include <dune/stuff/la/container/istl.hh>
#include <dune/stuff/common/parallel/helper.hh>

namespace Dune {
namespace GDT {
namespace Spaces {

template <class ViewImp,
          bool is_parallel =
              Dune::Stuff::UseParallelCommunication<typename ViewImp::Grid::CollectiveCommunication>::value>
struct CommunicationChooser
{
  typedef Dune::Stuff::SequentialCommunication Type;

  static Type* create(const ViewImp& /*gridView*/)
  {
    return new Type;
  }

  template <class SpaceBackend>
  static bool prepare(const SpaceBackend& /*space_backend*/, Type& /*communicator*/)
  {
    return false;
  }
}; // struct CommunicationChooser


#if HAVE_MPI


template <class ViewImp>
struct CommunicationChooser<ViewImp, true>
{
  typedef OwnerOverlapCopyCommunication<bigunsignedint<96>, int> Type;

  static Type* create(const ViewImp& gridView)
  {
    return new Type(gridView.comm());
  }

  template <class Space>
  static bool prepare(const Space& space, Type& communicator)
  {
    Stuff::LA::IstlRowMajorSparseMatrix<typename Space::RangeFieldType> matrix;
    PDELab::istl::ParallelHelper<typename Space::BackendType>(space.backend(), 0)
        .createIndexSetAndProjectForAMG(matrix.backend(), communicator);
    return true;
  }
}; // struct CommunicationChooser< ..., true >


#endif // HAVE_MPI

} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PARALLEL_HH
