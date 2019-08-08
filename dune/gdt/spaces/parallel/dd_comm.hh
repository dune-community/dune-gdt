// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   René Fritze     (2019)

#ifndef DUNE_GDT_SPACES_PARALLEL_DD_COMM_HH
#define DUNE_GDT_SPACES_PARALLEL_DD_COMM_HH


#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/istl/owneroverlapcopy.hh>

#include <dune/xt/common/math.hh>

#include <dune/xt/la/container/istl.hh>

#include <dune/gdt/tools/dd_helper.hh>

#include "dd_handles.hh"

namespace Dune {
namespace GDT {


template <class DdGridType> //, size_t r, size_t rD, class R>
class DomainDecompositionParallelHelper
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
  //      using SpaceType = SpaceInterface<GV, r, rD, R>;
  using GV = typename DdGridType::MacroGridViewType;

public:
  using DofCommunicatorType = typename DofCommunicationChooser<GV, true>::Type;

  DomainDecompositionParallelHelper(const DdGridType& dd_grid, std::string local_space_type, int verbose = 1)
    : dd_grid_(dd_grid)
    , local_space_type_(local_space_type)
    , rank_(dd_grid.macro_grid_view().comm().rank())
    // subdomain count ist nicht global summiert?!
    , rank_vector_(dd_grid.num_subdomains(), rank_)
    , ghosts_(dd_grid.num_subdomains(), false)
    , verbose_(verbose)
  {
    auto view = dd_grid_.macro_grid_view();

    // not optimal
    _interiorBorder_all_interface = InteriorBorder_All_Interface;
    _all_all_interface = All_All_Interface;

    if (view.comm().size() > 1) {
      // find out about ghosts
      GDT::DD::GhostDataHandle<DdGridType, GhostVector> gdh(dd_grid_, ghosts_, false);
      view.communicate(gdh, _interiorBorder_all_interface, Dune::ForwardCommunication);

      // create disjoint DOF partitioning
      //            SpaceTypeDataHandle<SpaceType,RankVector,DisjointPartitioningGatherScatter<RankIndex> >
      //  ibdh(dd_grid_,rank_vector_,DisjointPartitioningGatherScatter<RankIndex>(rank_));
      GDT::DD::DisjointPartitioningDataHandle<DdGridType, RankVector> pdh(dd_grid_, rank_vector_);
      view.communicate(pdh, _interiorBorder_all_interface, Dune::ForwardCommunication);
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
  const DdGridType& dd_grid_;
  const std::string local_space_type_;
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

template <class DdGridType>
void DomainDecompositionParallelHelper<DdGridType>::setup_parallel_indexset(DofCommunicatorType& dof_communicator)
{

  // ********************************************************************************
  // In the following, the code will always assume that all DOFs stored in a single
  // block of the BCRSMatrix are attached to the same entity and can be handled
  // identically. For that reason, the code often restricts itself to inspecting the
  // first entry of the blocks in the diverse BlockVectors.
  // ********************************************************************************

  const auto& macro_view = dd_grid_.macro_grid_view();
  using LocalSpaceType = std::unique_ptr<GDT::SpaceInterface<typename DdGridType::MicroGridViewType>>;
  using MacroEntity = typename DdGridType::MacroEntityType;
  using MacroId = typename DdGridType::MacroGridType::GlobalIdSetType::IdType;
  std::map<MacroId, LocalSpaceType> local_spaces;
  const auto& macro_id_set = macro_view.grid().globalIdSet();
  const auto& subdomain_id_to_entity_idx = dd_grid_.subdomain_id_to_idx_map();
  std::map<MacroId, size_t> local_sizes;
  size_t sum_local_sizes = 0;
  for (auto&& macro_element : elements(dd_grid_.macro_grid_view())) {
    typename DdGridType::MicroGridViewType view{dd_grid_.local_grid(macro_element).leaf_view()};
    LocalSpaceType lsp{make_subdomain_space(view, local_space_type_)};
    const auto id{macro_id_set.id(macro_element)};
    local_spaces.insert({id, std::move(lsp)});
    const size_t size{local_spaces.at(id)->mapper().size()};
    sum_local_sizes += size;
    local_sizes.insert({id, size});
  }
  // WAS space.mapper().size()
  const auto vector_size = sum_local_sizes;

  // Do we need to communicate at all?
  auto& comm = dd_grid_.global_comm();
  const bool need_communication = comm.size() > 1;

  // First find out which dofs we share with other processors
  BoolVector sharedDOF(vector_size, false);

  if (need_communication) {
    //    GDT::DD::SharedDOFDataHandle<DdGridType, BoolVector> data_handle(dd_grid_, sharedDOF, false);
    //    macro_view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  }

  // Count shared dofs that we own
  using GlobalIndex = typename DofCommunicatorType::ParallelIndexSet::GlobalIndex;
  GlobalIndex count{0};

  for (size_t i = 0; i < sharedDOF.size(); ++i) {
    if (owned_for_amg(i) && sharedDOF[i])
      ++count;
  }

  // Communicate per-rank count of owned and shared DOFs to all processes.
  std::vector<GlobalIndex> counts(comm.size());
  comm.allgather(&count, 1, &(counts[0]));

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
  //  if (need_communication) {
  //    GDT::DD::MinDataHandle<DdGridType, GlobalIndexVector> data_handle(dd_grid_, scalarIndices);
  //    macro_view.communicate(data_handle, _interiorBorder_all_interface, Dune::ForwardCommunication);
  //  }

  // Setup the index set
  //  dof_communicator.indexSet().beginResize();
  //  for (size_t i = 0; i < scalarIndices.size(); ++i) {
  //    Dune::OwnerOverlapCopyAttributeSet::AttributeSet attr;
  //    if (scalarIndices[i] != std::numeric_limits<GlobalIndex>::max()) {
  //      // global index exist in index set
  //      if (owned_for_amg(i)) {
  //        // This dof is managed by us.
  //        attr = Dune::OwnerOverlapCopyAttributeSet::owner;
  //      } else {
  //        attr = Dune::OwnerOverlapCopyAttributeSet::copy;
  //      }
  //      dof_communicator.indexSet().add(scalarIndices[i],
  //                                      typename DofCommunicatorType::ParallelIndexSet::LocalIndex(i, attr));
  //    }
  //  }
  dof_communicator.indexSet().endResize();

  // Compute neighbors using communication
  std::set<int> neighbors;

  //  if (need_communication) {
  //    DD::SpaceNeighborDataHandle<DdGridType, int> data_handle(dd_grid_, rank_, neighbors);
  //    macro_view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  //  }

  dof_communicator.remoteIndices().setNeighbours(neighbors);
  dof_communicator.remoteIndices().template rebuild<false>();
}

#endif // HAVE_MPI

} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_SPACES_PARALLEL_DD_COMM_HH
