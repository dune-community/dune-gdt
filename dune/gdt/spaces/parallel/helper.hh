// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk (2017 - 2018)

#ifndef DUNE_GDT_SPACES_PARALLEL_HELPER_HH
#define DUNE_GDT_SPACES_PARALLEL_HELPER_HH

#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/istl/owneroverlapcopy.hh>

#include <dune/xt/la/container/istl.hh>

#include "datahandles.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r, size_t rD, class R>
class GenericParallelHelper
{
  //! Type for storing rank values.
  using RankIndex = int;

  template <class S>
  using CommunicatedVector = XT::LA::IstlDenseVector<S>;
  //! Type used to store owner rank values of all DOFs.
  using RankVector = CommunicatedVector<RankIndex>;
  //! using bool as scalar type precludes using commondense
  using BoolVector = CommunicatedVector<uint_fast8_t>;
  //! Type used to store ghost flag of all DOFs,
  using GhostVector = BoolVector;
  using SpaceType = SpaceInterface<GV, r, rD, R>;
  using DofCommunicatorType = typename SpaceType::DofCommunicatorType;


public:
  GenericParallelHelper(const SpaceType& space, int verbose = 1)
    : space_(space)
    , rank_(space.grid_view().comm().rank())
    , rank_vector_(space.mapper().size(), rank_)
    , ghosts_(space.mapper().size(), false)
    , verbose_(verbose)
  {
    auto view = space.grid_view();

    // not optimal
    _interiorBorder_all_interface = InteriorBorder_All_Interface;
    _all_all_interface = All_All_Interface;

    if (view.comm().size() > 1) {
      // find out about ghosts
      GDT::GhostDataHandle<GV, r, rD, R, GhostVector> gdh(space, ghosts_, false);
      space.grid_view().communicate(gdh, _interiorBorder_all_interface, Dune::ForwardCommunication);

      // create disjoint DOF partitioning
      //            SpaceTypeDataHandle<SpaceType,RankVector,DisjointPartitioningGatherScatter<RankIndex> >
      //  ibdh(space_,rank_vector_,DisjointPartitioningGatherScatter<RankIndex>(rank_));
      GDT::DisjointPartitioningDataHandle<GV, r, rD, R, RankVector> pdh(space, rank_vector_);
      space.grid_view().communicate(pdh, _interiorBorder_all_interface, Dune::ForwardCommunication);
    }
  }

public:
  //! Returns the MPI rank of this process.
  RankIndex rank() const
  {
    return rank_;
  }

#if HAVE_MPI

  /* \brief Sets up the parallel communication information for AMG.
   *
   * \warning  This function silenty assumes that the matrix only has a single level
   *           of blocking and will not work correctly otherwise. Also note that AMG
   *           will only work correctly for P1 discretisations.
   *
   * \param dof_communicator  The parallel information object providing index set, interfaces and
   *           communicators.
   */
  void setup_parallel_indexset(DofCommunicatorType& dof_communicator);

private:
  // Checks whether a matrix block is owned by the current process. Used for the AMG
  // construction and thus assumes a single level of blocking and blocks with ownership
  // restricted to a single DOF.
  bool owned_for_amg(std::size_t i) const
  {
    return rank_vector_[i] == rank_;
  }

#endif // HAVE_MPI

private:
  const SpaceType& space_;
  const RankIndex rank_;
  RankVector rank_vector_; // vector to identify unique decomposition
  GhostVector ghosts_; // vector to identify ghost dofs
  int verbose_;

  //! The actual communication interface used when algorithm requires InteriorBorder_All_Interface.
  InterfaceType _interiorBorder_all_interface;

  //! The actual communication interface used when algorithm requires All_All_Interface.
  InterfaceType _all_all_interface;
};

#if HAVE_MPI

template <class GV, size_t r, size_t rD, class R>
void GenericParallelHelper<GV, r, rD, R>::setup_parallel_indexset(DofCommunicatorType& dof_communicator)
{

  // ********************************************************************************
  // In the following, the code will always assume that all DOFs stored in a single
  // block of the BCRSMatrix are attached to the same entity and can be handled
  // identically. For that reason, the code often restricts itself to inspecting the
  // first entry of the blocks in the diverse BlockVectors.
  // ********************************************************************************
  const auto& view = space_.grid_view();
  const auto vector_size = space_.mapper().size();

  // Do we need to communicate at all?
  const bool need_communication = view.comm().size() > 1;

  // First find out which dofs we share with other processors
  BoolVector sharedDOF(vector_size, false);

  if (need_communication) {
    GDT::SharedDOFDataHandle<GV, r, rD, R, BoolVector> data_handle(space_, sharedDOF, false);
    view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  }

  // Count shared dofs that we own
  using GlobalIndex = typename DofCommunicatorType::ParallelIndexSet::GlobalIndex;
  GlobalIndex count{0};

  for (size_t i = 0; i < sharedDOF.size(); ++i) {
    if (owned_for_amg(i) && sharedDOF[i])
      ++count;
  }

  // Communicate per-rank count of owned and shared DOFs to all processes.
  std::vector<GlobalIndex> counts(view.comm().size());
  view.comm().allgather(&count, 1, &(counts[0]));

  // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
  GlobalIndex start = std::accumulate(counts.begin(), counts.begin() + rank_, GlobalIndex(0));

  using GlobalIndexVector = CommunicatedVector<GlobalIndex>;
  GlobalIndexVector scalarIndices(vector_size, std::numeric_limits<GlobalIndex>::max());

  for (size_t i = 0; i < sharedDOF.size(); ++i) {
    if (owned_for_amg(i) && sharedDOF[i]) {
      scalarIndices[i] = start;
      ++start;
    }
  }

  // Publish global indices for the shared DOFS to other processors.
  if (need_communication) {
    GDT::MinDataHandle<GV, r, rD, R, GlobalIndexVector> data_handle(space_, scalarIndices);
    view.communicate(data_handle, _interiorBorder_all_interface, Dune::ForwardCommunication);
  }

  // Setup the index set
  dof_communicator.indexSet().beginResize();
  for (size_t i = 0; i < scalarIndices.size(); ++i) {
    Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
    if (scalarIndices[i] != std::numeric_limits<GlobalIndex>::max()) {
      // global index exist in index set
      if (owned_for_amg(i)) {
        // This dof is managed by us.
        attr = Dune::OwnerOverlapCopyAttributeSet::owner;
      } else {
        attr = Dune::OwnerOverlapCopyAttributeSet::copy;
      }
      dof_communicator.indexSet().add(scalarIndices[i],
                                      typename DofCommunicatorType::ParallelIndexSet::LocalIndex(i, attr));
    }
  }
  dof_communicator.indexSet().endResize();

  // Compute neighbors using communication
  std::set<int> neighbors;

  if (need_communication) {
    SpaceNeighborDataHandle<GV, r, rD, R, int> data_handle(space_, rank_, neighbors);
    view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  }

  dof_communicator.remoteIndices().setNeighbours(neighbors);
  dof_communicator.remoteIndices().template rebuild<false>();
}

#endif // HAVE_MPI

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PARALLEL_HELPER_HH
