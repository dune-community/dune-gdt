
#ifndef DUNE_GDT_PARALLEL_HELPER_HH
#define DUNE_GDT_PARALLEL_HELPER_HH

#include <limits>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/typetraits.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/io.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>

namespace Dune {
namespace GDT {
template <class SpaceTraits>
class PdelabDataHandleShroudTraits : public SpaceTraits
{
public:
  using size_type = std::size_t;
  using SizeType = size_type;
  using ContainerIndex = PDELab::SimpleContainerIndex<std::size_t>;
  using DOFIndex = PDELab::SimpleDOFIndex<std::size_t>;
};
template <class SpaceImp>
class PdelabDataHandleShroud //: public SpaceImp
{
  using ThisType = PdelabDataHandleShroud<SpaceImp>;

public:
  using Traits = PdelabDataHandleShroudTraits<typename SpaceImp::Traits>;
  // this way index types are extracted from our traits and no manufacturing of an entire Ordering struct is needed
  using NodeTag = TypeTree::LeafNodeTag;
  using Ordering = ThisType;

  using ContainerIndex = typename Traits::ContainerIndex;
  using DOFIndex = typename Traits::DOFIndex;
  using size_type = typename Traits::SizeType;

  using GridType = XT::Grid::extract_grid_t<typename SpaceImp::GridViewType>;
  using Layer = typename XT::Grid::Layer<GridType, XT::Grid::Layers::leaf, XT::Grid::Backends::view>;
  using RealGridView = typename Layer::type;
  using EntitySet = PDELab::PartitionViewEntitySet<RealGridView, Partitions::All>;

private:
  template <typename V, typename EntityIndex>
  void setup_dof_indices(V& v, size_type n, const EntityIndex& ei, std::integral_constant<bool, true>) const
  {
    v.resize(n);
    for (typename V::iterator it = v.begin(), endit = v.end(); it != endit; ++it) {
      it->treeIndex().clear();
      it->entityIndex() = ei;
    }
  }

  template <typename V, typename EntityIndex>
  void setup_dof_indices(V& v, size_type n, const EntityIndex& ei, std::integral_constant<bool, false>) const
  {
  }

public:
  PdelabDataHandleShroud(SpaceImp space)
    : space_(space)
    , real_view_(Layer::create(space.grid_view().grid()))
    , entity_set_(real_view_)
  {
  }

  const typename SpaceImp::GridViewType& gridView() const
  {
    return space_.grid_view();
  }

  const auto maxLocalSize() const
  {
    return space_.mapper().maxNumDofs();
  }

  //! returns true if data for this codim should be communicated
  bool dataHandleContains(int codim) const
  {
    //    return space().ordering().contains(codim);
    //    return codim == 0;
    return true;
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool dataHandleFixedSize(int /*codim*/) const
  {
    return true;
  }

  //! send the leaf ordering offsets as exported by the EntityIndexCache
  bool sendLeafSizes() const
  {
    return false;
  }

  //! return vector of global indices associated with the given entity
  template <typename Entity, typename ContainerIndex, typename DOFIndex, typename OffsetIterator, bool map_dof_indices>
  void dataHandleIndices(const Entity& e,
                         std::vector<ContainerIndex>& container_indices,
                         std::vector<DOFIndex>& dof_indices,
                         OffsetIterator oit,
                         std::integral_constant<bool, map_dof_indices> map_dof_indices_value) const
  {

    static_assert((std::is_same<ContainerIndex, typename Ordering::Traits::ContainerIndex>::value),
                  "dataHandleContainerIndices() called with invalid ContainerIndex type.");
    static const size_type entity_capacity = 1;
    using EntityIndex = std::array<size_type, 2>;
    EntityIndex ei;
    const size_type size{dataHandleSize(e)};
    ei[0] = GlobalGeometryTypeIndex::index(e.type());
    ei[1] = entity_set_.indexSet().index(e);
    //    DUNE_THROW(NotImplemented, "");
    //        PDELab::get_leaf_offsets_for_entity<EntityIndex,OffsetIterator> get_offsets(ei,oit);
    //    TypeTree::applyToTree(space().ordering(),get_offsets);
    //    OffsetIterator end_oit = oit + (TypeTree::TreeInfo<Ordering>::leafCount + 1);
    //
    //    // convert sizes to offsets - last entry contains total size
    //    std::partial_sum(oit,end_oit,oit);
    //    size_type size = *(oit + TypeTree::TreeInfo<Ordering>::leafCount);
    //
    //    container_indices.resize(size);
    //    // Clear index state
    //    for (typename std::vector<ContainerIndex>::iterator it = container_indices.begin(),
    //             endit = container_indices.end();
    //         it != endit;
    //         ++it)
    //      it->clear();
    //
    setup_dof_indices(dof_indices, size, ei, map_dof_indices_value);
    //

    const size_type geometry_type_index = ei[0];
    const size_type entity_index = ei[1];
    auto ci_out = container_indices.begin();
    for (size_type i = 0; i < size; ++i, ++ci_out) {
      //                    ci_out->push_back(i);
      *ci_out += entity_index * size;
    }
  }

  template <typename Entity>
  size_type dataHandleSize(const Entity& e) const
  {
    return e.subEntities(2);
  }

  const EntitySet& entitySet() const
  {
    return entity_set_;
  }

private:
  SpaceImp space_;
  const RealGridView real_view_;
  const EntitySet entity_set_;
};

template <class S>
struct is_space<PdelabDataHandleShroud<S>, false> : public std::true_type
{
};


//========================================================
// A parallel helper class providing a nonoverlapping
// decomposition of all degrees of freedom
//========================================================

template <typename SpaceType>
class GenericParallelHelper
{

  //! Type for storing rank values.
  using RankIndex = int;

  //! Type used to store owner rank values of all DOFs.
  using RankVector = XT::LA::IstlDenseVector<RankIndex>;
  //! Type used to store ghost flag of all DOFs.
  using GhostVector = XT::LA::IstlDenseVector<bool>;
  //! this adds necesasary types for DataHandles
  using ShroudedSpaceType = PdelabDataHandleShroud<SpaceType>;
  using IndexCache = XT::LA::XTEntityIndexCache<ShroudedSpaceType>;

public:
  GenericParallelHelper(const SpaceType& space, int verbose = 1)
    : space_(space)
    , rank_(space.grid_view().comm().rank())
    , rank_vector_(space.grid_view().grid().size(0), rank_)
    , ghosts_(space.grid_view().grid().size(0), false)
    , verbose_(verbose)
  {
    using GridType = std::decay_t<decltype(space.grid_view().grid())>;
    auto view =
        XT::Grid::Layer<GridType, XT::Grid::Layers::leaf, XT::Grid::Backends::view>::create(space.grid_view().grid());
    using GridViewType = decltype(view);
    PDELab::impl::EntitySet<GridViewType> entity_set(view);
    // Let's try to be clever and reduce the communication overhead by picking the smallest
    // possible communication interface depending on the overlap structure of the SpaceType.
    // FIXME: Switch to simple comparison as soon as dune-grid:1b3e83ec0 is reliably available!
    if (entity_set.partitions().value == Partitions::interiorBorder.value) {
      // The SpaceType only spans the interior and border partitions, so we can skip sending or
      // receiving anything else.
      _interiorBorder_all_interface = InteriorBorder_InteriorBorder_Interface;
      _all_all_interface = InteriorBorder_InteriorBorder_Interface;
    } else {
      // In general, we have to transmit more.
      _interiorBorder_all_interface = InteriorBorder_All_Interface;
      _all_all_interface = All_All_Interface;
    }

    //! too broad, see aboce
    _interiorBorder_all_interface = InteriorBorder_All_Interface;
    _all_all_interface = All_All_Interface;

    if (view.comm().size() > 1) {

      ShroudedSpaceType shr(space_);
      // find out about ghosts
      DXTC_LOG_DEBUG << "local view size " << space.grid_view().grid().size(0) << std::endl;
      PDELab::GhostDataHandle<ShroudedSpaceType, GhostVector, IndexCache> gdh(shr, ghosts_, false);
      space.grid_view().communicate(gdh, _interiorBorder_all_interface, Dune::ForwardCommunication);
      DXTC_LOG_DEBUG << "GHOSTS " << ghosts_ << std::endl;

      // create disjoint DOF partitioning
      //            SpaceTypeDataHandle<SpaceType,RankVector,DisjointPartitioningGatherScatter<RankIndex> >
      //  ibdh(space_,rank_vector_,DisjointPartitioningGatherScatter<RankIndex>(rank_));
      PDELab::DisjointPartitioningDataHandle<ShroudedSpaceType, RankVector, IndexCache> pdh(shr, rank_vector_);
      space.grid_view().communicate(pdh, _interiorBorder_all_interface, Dune::ForwardCommunication);
      //      DXTC_LOG_DEBUG << "RANKS " << rank_vector_ << std::endl;
    }
  }


public:
  //! Tests whether the given index is owned by this process.
  template <class ContainerIndex>
  bool owned(const ContainerIndex& i) const
  {
    return rank_vector_[i] == rank_;
  }

  //! Tests whether the given index belongs to a ghost DOF.
  template <class ContainerIndex>
  bool isGhost(const ContainerIndex& i) const
  {
    return ghosts_[i];
  }

public:
  //! Returns the MPI rank of this process.
  RankIndex rank() const
  {
    return rank_;
  }

#if HAVE_MPI

  //! Makes the matrix consistent and creates the parallel information for AMG.
  /**
   * This function accomplishes two things:
   *
   * 1. Makes the matrix consistent w.r.t. to the disjoint partitioning of the DOF space,
   *    i.e. aggregates matrix entries for border entries from neighboring ranks.
   *
   * 2. Sets up the parallel communication information for AMG.
   *
   * \warning  This function silenty assumes that the matrix only has a single level
   *           of blocking and will not work correctly otherwise. Also note that AMG
   *           will only work correctly for P1 discretisations.
   *
   * \param m  The PDELab matrix container.
   * \param c  The parallel information object providing index set, interfaces and
   *           communicators.
   */
  template <typename MatrixType, typename Comm>
  void createIndexSetAndProjectForAMG(MatrixType& m, Comm& c);

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

template <typename SpaceType>
template <typename M, typename C>
void GenericParallelHelper<SpaceType>::createIndexSetAndProjectForAMG(M& m, C& c)
{
  const bool is_bcrs_matrix =
      std::is_same<typename PDELab::ISTL::tags::container<PDELab::Backend::Native<M>>::type::base_tag,
                   PDELab::ISTL::tags::bcrs_matrix>::value;

  const bool block_type_is_field_matrix = std::is_same<
      typename PDELab::ISTL::tags::container<typename PDELab::Backend::Native<M>::block_type>::type::base_tag,
      PDELab::ISTL::tags::field_matrix>::value;

  // We assume M to be a BCRSMatrix in the following, so better check for that
  static_assert(is_bcrs_matrix && block_type_is_field_matrix, "matrix structure not compatible with AMG");

  // ********************************************************************************
  // In the following, the code will always assume that all DOFs stored in a single
  // block of the BCRSMatrix are attached to the same entity and can be handled
  // identically. For that reason, the code often restricts itself to inspecting the
  // first entry of the blocks in the diverse BlockVectors.
  // ********************************************************************************

  const auto& view = space_.grid_view();
  const auto vector_size = space_.grid_view().grid().size(0);
  ShroudedSpaceType shr(space_);

  // Do we need to communicate at all?
  const bool need_communication = view.comm().size() > 1;

  // First find out which dofs we share with other processors
  using BoolVector = XT::LA::IstlDenseVector<bool>;
  BoolVector sharedDOF(vector_size, false);

  if (need_communication) {
    PDELab::SharedDOFDataHandle<ShroudedSpaceType, BoolVector, IndexCache> data_handle(shr, sharedDOF, false);
    view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  }

  // Count shared dofs that we own
  typedef typename C::ParallelIndexSet::GlobalIndex GlobalIndex;
  GlobalIndex count = GlobalIndex(0);

  for (size_t i = 0; i < sharedDOF.size(); ++i)
    if (owned_for_amg(i) && sharedDOF[i])
      ++count;

  DXTC_LOG_DEBUG << view.comm().rank() << ": shared block count is " << count << std::endl;

  // Communicate per-rank count of owned and shared DOFs to all processes.
  std::vector<GlobalIndex> counts(view.comm().size());
  view.comm().allgather(&count, 1, &(counts[0]));

  // Compute start index start_p = \sum_{i=0}^{i<p} counts_i
  GlobalIndex start = std::accumulate(counts.begin(), counts.begin() + rank_, GlobalIndex(0));

  using GIVector = XT::LA::IstlDenseVector<GlobalIndex>;
  GIVector scalarIndices(vector_size, std::numeric_limits<GlobalIndex>::max());

  for (size_t i = 0; i < sharedDOF.size(); ++i)
    if (owned_for_amg(i) && sharedDOF[i]) {
      scalarIndices[i] = start;
      ++start;
    }

  // Publish global indices for the shared DOFS to other processors.
  if (need_communication) {
    PDELab::MinDataHandle<ShroudedSpaceType, GIVector, IndexCache> data_handle(shr, scalarIndices);
    view.communicate(data_handle, _interiorBorder_all_interface, Dune::ForwardCommunication);
  }

  // Setup the index set
  c.indexSet().beginResize();
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
      c.indexSet().add(scalarIndices[i], typename C::ParallelIndexSet::LocalIndex(i, attr));
    }
  }
  c.indexSet().endResize();

  // Compute neighbors using communication
  std::set<int> neighbors;

  if (need_communication) {
    PDELab::GFSNeighborDataHandle<ShroudedSpaceType, int> data_handle(shr, rank_, neighbors);
    view.communicate(data_handle, _all_all_interface, Dune::ForwardCommunication);
  }

  c.remoteIndices().setNeighbours(neighbors);
  c.remoteIndices().template rebuild<false>();
}

#endif // HAVE_MPI

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PARALLEL_HELPER_HH
