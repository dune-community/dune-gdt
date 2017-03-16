// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2014, 2016)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_PARALLEL_HH
#define DUNE_GDT_SPACES_PARALLEL_HH

#include <dune/common/parallel/communicator.hh>

#include <dune/istl/owneroverlapcopy.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/backend/istl.hh>
#endif

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/common/parallel/helper.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/spaces/parallel_helper.hh>

namespace Dune {
namespace GDT {


template <class ViewImp,
          bool is_parallel = Dune::XT::UseParallelCommunication<typename ViewImp::Grid::CollectiveCommunication>::value>
struct CommunicationChooser
{
  typedef Dune::XT::SequentialCommunication Type;

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
  using type = Type;

  static Type* create(const ViewImp& gridView)
  {
    return new Type(gridView.comm());
  }

  template <class Space>
  static typename std::enable_if<Space::backend_type == Dune::GDT::ChooseSpaceBackend::pdelab, bool>::type
  prepare(const Space& space, Type& communicator)
  {
#if HAVE_DUNE_PDELAB
    XT::LA::IstlRowMajorSparseMatrix<typename Space::RangeFieldType> matrix;
    PDELab::istl::ParallelHelper<typename Space::BackendType>(space.backend(), 0)
        .createIndexSetAndProjectForAMG(matrix.backend(), communicator);
#endif // HAVE_DUNE_PDELAB
    return true;
  } // ... prepare(...)

  template <class Space>
  static typename std::enable_if<Space::backend_type == Dune::GDT::ChooseSpaceBackend::gdt, bool>::type
  prepare(const Space& /*space*/, Type& communicator)
  {
    GDT::GenericParallelHelper<Space>(space, 1).createIndexSetAndProjectForAMG(communicator);
    return true;
  } // ... prepare(...)

  template <class Space>
  static typename std::enable_if<Space::backend_type == Dune::GDT::ChooseSpaceBackend::fem, bool>::type
  prepare(const Space& /*space*/, Type& communicator)
  {
    communicator.remoteIndices().template rebuild<true>();
    return true;
  } // ... prepare(...)

}; // struct CommunicationChooser< ..., true >


#endif // HAVE_MPI

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PARALLEL_HH
