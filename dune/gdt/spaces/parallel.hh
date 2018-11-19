// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2018)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_PARALLEL_HH
#define DUNE_GDT_SPACES_PARALLEL_HH

#include <dune/common/parallel/communicator.hh>

#include <dune/istl/owneroverlapcopy.hh>

#include <dune/xt/la/container/istl.hh>
#include <dune/xt/common/parallel/helper.hh>
//#include <dune/xt/grid/layers.hh>

#include <dune/gdt/spaces/parallel_helper.hh>

namespace Dune {
namespace GDT {


template <class ViewImp,
          bool is_parallel = Dune::XT::UseParallelCommunication<
              typename XT::Grid::extract_grid<ViewImp>::type::CollectiveCommunication>::value>
struct DofCommunicationChooser
{
  typedef Dune::XT::SequentialCommunication Type;
  static_assert(Dune::AlwaysFalse<ViewImp>::value, "");

  static Type* create(const ViewImp& /*gridView*/)
  {
    static_assert(Dune::AlwaysFalse<ViewImp>::value, "");
    return new Type;
  }

  template <class SpaceBackend>
  static bool prepare(const SpaceBackend& /*space_backend*/, Type& /*communicator*/)
  {
    static_assert(Dune::AlwaysFalse<SpaceBackend>::value, "");
    return false;
  }
}; // struct DofCommunicationChooser


#if HAVE_MPI


template <class ViewImp>
struct DofCommunicationChooser<ViewImp, true>
{
private:
  // this is necessary because alugrid's id is not integral
  using RealGlobalId = typename XT::Grid::extract_grid_t<ViewImp>::GlobalIdSet::IdType;
  using RealLocalId = typename XT::Grid::extract_grid_t<ViewImp>::LocalIdSet::IdType;

public:
  using GlobalId = std::conditional_t<std::is_arithmetic<RealGlobalId>::value, RealGlobalId, uint64_t>;
  using LocalId = std::conditional_t<std::is_arithmetic<RealLocalId>::value, RealLocalId, int>;
  using Type = OwnerOverlapCopyCommunication<GlobalId, LocalId>;
  using type = Type;

  static Type* create(const ViewImp& gridView)
  {
    return new Type(gridView.comm(), SolverCategory::overlapping);
  }

  template <class Space>
  static bool prepare(const Space& space, Type& communicator)
  {
    DXTC_LOG_DEBUG << "Preparing dof_comm ...\n";
    GDT::GenericParallelHelper<Space>(space, 1).setup_parallel_indexset(communicator);
    DXTC_LOG_DEBUG << "Preparing dof_comm done\n";
    return true;
  } // ... prepare(...)

}; // struct DofCommunicationChooser< ..., true >


#endif // HAVE_MPI

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PARALLEL_HH
