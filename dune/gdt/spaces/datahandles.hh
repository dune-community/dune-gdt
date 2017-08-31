#ifndef DUNE_GDT_DATAHANDLES_HH
#define DUNE_GDT_DATAHANDLES_HH

#include <dune/gdt/spaces/localview.hh>

namespace Dune {
namespace GDT {

//! Communication descriptor for sending one item of type E per DOF.
//! used in Mindatahandle
template <typename E>
struct DOFDataCommunicationDescriptor
{
  typedef E DataType;

  // Wrap the grid's communication buffer to enable sending leaf ordering sizes along with the data
  static const bool wrap_buffer = true;

  // export original data type to fix up size information forwarded to standard gather / scatter functors
  typedef E OriginalDataType;

  template <typename SpaceType>
  bool contains(const SpaceType& space, int dim, int codim) const
  {
    return SpaceType::associates_data_with(codim);
  }

  template <typename SpaceType>
  bool fixedSize(const SpaceType& /*space*/, int /*dim*/, int /*codim*/) const
  {
    // we currently have no adaptive spaces
    return true;
  }

  template <typename SpaceType, typename Entity>
  std::size_t size(const SpaceType& space, const Entity& e) const
  {
    return space.grid_layer().indexSet().contains(e) ? space.mapper().numDofs(e) : 0u;
  }
};

//! Communication descriptor for sending count items of type E per entity with attached DOFs.
//! used in  SharedDOFDataHandle, DisjointPartitioningDataHandle, ghostdata
template <typename E>
struct EntityDataCommunicationDescriptor
{
  typedef E DataType;
  // Data is per entity, so we don't need to send leaf ordering size and thus can avoid wrapping the
  // grid's communication buffer
  static const bool wrap_buffer = false;

  template <typename SpaceType>
  bool contains(const SpaceType& /*space*/, int /*dim*/, int codim) const
  {
    return SpaceType::associates_data_with(codim);
  }

  template <typename SpaceType>
  bool fixedSize(const SpaceType& /*space*/, int /*dim*/, int /*codim*/) const
  {
    return true;
  }

  template <typename SpaceType, typename Entity>
  std::size_t size(const SpaceType& space, const Entity& e) const
  {
    return (SpaceType::associates_data_with(Entity::codimension) && space.grid_layer().indexSet().contains(e)) ? _count
                                                                                                               : 0;
  }

  //! remove default value, force handles to use actual space provided values
  explicit EntityDataCommunicationDescriptor(std::size_t count = 4)
    : _count(count)
  {
  }

private:
  const std::size_t _count;
};


template <typename SpaceType,
          typename VectorType,
          typename GatherScatter,
          typename CommunicationDescriptor = DOFDataCommunicationDescriptor<typename VectorType::ScalarType>>
class SpaceDataHandle
    : public Dune::CommDataHandleIF<SpaceDataHandle<SpaceType, VectorType, GatherScatter, CommunicationDescriptor>,
                                    typename CommunicationDescriptor::DataType>
{

public:
  using DataType = typename CommunicationDescriptor::DataType;
  using Codim0EntityType = XT::Grid::extract_entity_t<typename SpaceType::GridViewType>;

  SpaceDataHandle(const SpaceType& space,
                  VectorType& v,
                  GatherScatter gather_scatter = GatherScatter(),
                  CommunicationDescriptor communication_descriptor = CommunicationDescriptor())
    : space_(space)
    , _gather_scatter(gather_scatter)
    , _communication_descriptor(communication_descriptor)
    , local_view_(v, space, _communication_descriptor)
  {
  }

  //! returns true if data for this codim should be communicated
  bool contains(int dim, int codim) const
  {
    return _communication_descriptor.contains(space_, dim, codim);
  }

  //!  \brief returns true if size per entity of given dim and codim is a constant
  bool fixedsize(int dim, int codim) const
  {
    return _communication_descriptor.fixedSize(space_, dim, codim);
  }

  /*!  \brief how many objects of type DataType have to be sent for a given entity

    Note: Only the sender side needs to know this size.
  */
  template <class EntityType>
  size_t size(const EntityType& entity) const
  {
    return _communication_descriptor.size(space_, entity);
  }

  /** pack data from user to message buffer
   *  the parent type prescribes a CONST gather method,
   *    because obviously receiving data will never change our local state
   *  this needs to be defined for all codims, but only codims for which
   *  contains is true will ever be called at runtime
   */
  template <typename MessageBuffer, class EntityType>
  void gather(MessageBuffer& buffer, const EntityType& entity) const
  {
    local_view_.bind(entity);
    if (_gather_scatter.gather(buffer, entity, local_view_))
      local_view_.commit();
  }

  /** unpack data from message buffer to user
    n is the number of objects sent by the sender
   *  this needs to be defined for all codims, but only codims for which
   *  contains is true will ever be called at runtime
  */
  template <typename MessageBuffer, class EntityType>
  void scatter(MessageBuffer& buff, const EntityType& e, size_t n)
  {
    local_view_.bind(e);
    if (_gather_scatter.scatter(buff, n, e, local_view_))
      local_view_.commit();
  }


private:
  const SpaceType& space_;
  //! has to be mutable because gather is const
  mutable GatherScatter _gather_scatter;
  CommunicationDescriptor _communication_descriptor;
  mutable LocalView<typename VectorType::Traits, typename VectorType::ScalarType, SpaceType, CommunicationDescriptor>
      local_view_;
};


template <typename GatherScatter>
class DataGatherScatter
{

public:
  typedef std::size_t size_t;

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool gather(MessageBuffer& buff, const Entity& e, const LocalView& local_view) const
  {
    for (std::size_t i = 0; i < local_view.size(); ++i)
      _gather_scatter.gather(buff, local_view[i]);
    return false;
  }

  // default scatter - requires function space structure to be identical on sender and receiver side
  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool scatter(MessageBuffer& buff, size_t n, const Entity& e, LocalView& local_view) const
  {
    if (partitions_.contains(e.partitionType())) {
      if (local_view.size() != n)
        DUNE_THROW(Exception,
                   "size mismatch in GridFunctionSpace data handle, have " << local_view.size() << "DOFs, but received "
                                                                           << n);

      for (std::size_t i = 0; i < local_view.size(); ++i)
        _gather_scatter.scatter(buff, local_view[i]);
      return true;
    } else {
      if (local_view.size() != 0)
        DUNE_THROW(Exception,
                   "expected no DOFs in partition '" << e.partitionType() << "', but have " << local_view.size());

      for (std::size_t i = 0; i < local_view.size(); ++i) {
        typename LocalView::value_type dummy;
        buff.read(dummy);
      }
      return false;
    }
  }

  DataGatherScatter(GatherScatter gather_scatter = GatherScatter())
    : _gather_scatter(gather_scatter)
  {
  }

private:
  GatherScatter _gather_scatter;
  // this is currently the default since we have no non-overlapping parallel code
  static const Partitions::All partitions_;
};


template <typename GatherScatter>
class DataEntityGatherScatter
{

public:
  typedef std::size_t size_t;

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool gather(MessageBuffer& buff, const Entity& e, const LocalView& local_view) const
  {
    for (std::size_t i = 0; i < local_view.size(); ++i)
      _gather_scatter.gather(buff, e, local_view[i]);
    return false;
  }

  // see documentation in DataGatherScatter for further info on the scatter() implementations
  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool scatter(MessageBuffer& buff, size_t n, const Entity& e, LocalView& local_view) const
  {
    if (local_view.cache().gridFunctionSpace().entitySet().partitions().contains(e.partitionType())) {
      if (local_view.size() != n)
        DUNE_THROW(Exception,
                   "size mismatch in GridFunctionSpace data handle, have " << local_view.size() << "DOFs, but received "
                                                                           << n);

      for (std::size_t i = 0; i < local_view.size(); ++i)
        _gather_scatter.scatter(buff, e, local_view[i]);
      return true;
    } else {
      if (local_view.size() != 0)
        DUNE_THROW(Exception,
                   "expected no DOFs in partition '" << e.partitionType() << "', but have " << local_view.size());

      for (std::size_t i = 0; i < local_view.size(); ++i) {
        typename LocalView::ElementType dummy;
        buff.read(dummy);
      }
      return false;
    }
  }

  DataEntityGatherScatter(GatherScatter gather_scatter = GatherScatter())
    : _gather_scatter(gather_scatter)
  {
  }

private:
  GatherScatter _gather_scatter;
};

class MinGatherScatter
{
public:
  template <class MessageBuffer, class DataType>
  void gather(MessageBuffer& buff, DataType& data) const
  {
    buff.write(data);
  }

  template <class MessageBuffer, class DataType>
  void scatter(MessageBuffer& buff, DataType& data) const
  {
    DataType x;
    buff.read(x);
    data = std::min(data, x);
  }
};

template <class SpaceType, class VectorType>
class MinDataHandle : public SpaceDataHandle<SpaceType,
                                             VectorType,
                                             DataGatherScatter<MinGatherScatter>,
                                             DOFDataCommunicationDescriptor<typename VectorType::ScalarType>>
{
  typedef SpaceDataHandle<SpaceType,
                          VectorType,
                          DataGatherScatter<MinGatherScatter>,
                          DOFDataCommunicationDescriptor<typename VectorType::ScalarType>>
      BaseType;

public:
  MinDataHandle(const SpaceType& space_, VectorType& v_)
    : BaseType(space_, v_)
  {
  }
};

//! GatherScatter functor for marking ghost DOFs.
/**
 * This data handle will mark all ghost DOFs (more precisely, all DOFs associated
 * with entities not part of either the interior or the border partition).
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::InteriorBorder_All_Interface.
 */
class GhostGatherScatter
{
public:
  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool gather(MessageBuffer& buff, const Entity& e, LocalView& local_view) const
  {
    // Figure out where we are...
    const bool ghost = e.partitionType() != Dune::InteriorEntity && e.partitionType() != Dune::BorderEntity;

    // ... and send something (doesn't really matter what, we'll throw it away on the receiving side).
    buff.write(ghost);

    return false;
  }

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool scatter(MessageBuffer& buff, std::size_t n, const Entity& e, LocalView& local_view) const
  {
    // Figure out where we are - we have to do this again on the receiving side due to the asymmetric
    // communication interface!
    const bool ghost = e.partitionType() != Dune::InteriorEntity && e.partitionType() != Dune::BorderEntity;

    // drain buffer
    bool dummy;
    buff.read(dummy);

    for (std::size_t i = 0; i < local_view.size(); ++i)
      local_view[i] = ghost;

    return true;
  }
};

//! Data handle for marking ghost DOFs.
/**
 * This data handle will mark all ghost DOFs (more precisely, all DOFs associated
 * with entities not part of either the interior or the border partition).
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::InteriorBorder_All_Interface.
 */
template <class SpaceType, class VectorType>
class GhostDataHandle
    : public SpaceDataHandle<SpaceType, VectorType, GhostGatherScatter, EntityDataCommunicationDescriptor<bool>>
{
  typedef SpaceDataHandle<SpaceType, VectorType, GhostGatherScatter, EntityDataCommunicationDescriptor<bool>> BaseType;

  //  static_assert((std::is_same<typename VectorType::ScalarType, bool>::value),
  //                "GhostDataHandle expects a vector of bool values");

public:
  //! Creates a new GhostDataHandle.
  /**
   * Creates a new GhostDataHandle and by default initializes the result vector
   * with the correct value of false. If you have already done that externally,
   * you can skip the initialization.
   *
   * \param space_         The GridFunctionSpace to operate on.
   * \param v_           The result vector.
   * \param init_vector  Flag to control whether the result vector will be initialized.
   */
  GhostDataHandle(const SpaceType& space_, VectorType& v_, bool init_vector = true)
    : BaseType(space_, v_)
  {
    if (init_vector)
      v_.set_all(false);
  }
};


//! GatherScatter functor for creating a disjoint DOF partitioning.
/**
 * This functor will associate each DOF with a unique rank, creating a nonoverlapping partitioning
 * of the unknowns. The rank for a DOF is chosen by finding the lowest rank on which the associated
 * grid entity belongs to either the interior or the border partition.
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::InteriorBorder_All_Interface and the result vector must be initialized with the MPI rank value.
 */
template <typename RankIndex>
class DisjointPartitioningGatherScatter
{

public:
  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool gather(MessageBuffer& buff, const Entity& e, LocalView& local_view) const
  {
    // We only gather from interior and border entities, so we can throw in our ownership
    // claim without any further checks.
    buff.write(_rank);

    return true;
  }

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool scatter(MessageBuffer& buff, std::size_t n, const Entity& e, LocalView& local_view) const
  {
    // Value used for DOFs with currently unknown rank.
    const RankIndex unknown_rank = std::numeric_limits<RankIndex>::max();

    // We can only own this DOF if it is either on the interior or border partition.
    const bool is_interior_or_border =
        (e.partitionType() == Dune::InteriorEntity || e.partitionType() == Dune::BorderEntity);

    // Receive data.
    RankIndex received_rank;
    buff.read(received_rank);

    for (std::size_t i = 0; i < local_view.size(); ++i) {
      // Get the currently stored owner rank for this DOF.
      RankIndex current_rank = local_view[i];

      // We only gather from interior and border entities, so we need to make sure
      // we relinquish any ownership claims on overlap and ghost entities on the
      // receiving side. We also need to make sure not to overwrite any data already
      // received, so we only blank the rank value if the currently stored value is
      // equal to our own rank.
      if (!is_interior_or_border && current_rank == _rank)
        current_rank = unknown_rank;

      // Assign DOFs to minimum rank value.
      local_view[i] = std::min(current_rank, received_rank);
    }
    return true;
  }

  //! Create a DisjointPartitioningGatherScatter object.
  /**
   * \param rank  The MPI rank of the current process.
   */
  DisjointPartitioningGatherScatter(RankIndex rank)
    : _rank(rank)
  {
  }

private:
  const RankIndex _rank;
};

//! GatherScatter data handle for creating a disjoint DOF partitioning.
/**
 * This data handle will associate each DOF with a unique rank, creating a nonoverlapping partitioning
 * of the unknowns. The rank for a DOF is chosen by finding the lowest rank on which the associated
 * grid entity belongs to either the interior or the border partition.
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::InteriorBorder_All_Interface and the result vector must be initialized with the MPI rank value.
 */
template <class SpaceType, class VectorType>
class DisjointPartitioningDataHandle
    : public SpaceDataHandle<SpaceType,
                             VectorType,
                             DisjointPartitioningGatherScatter<typename VectorType::ScalarType>,
                             EntityDataCommunicationDescriptor<typename VectorType::ScalarType>>
{
  typedef SpaceDataHandle<SpaceType,
                          VectorType,
                          DisjointPartitioningGatherScatter<typename VectorType::ScalarType>,
                          EntityDataCommunicationDescriptor<typename VectorType::ScalarType>>
      BaseType;

public:
  //! Creates a new DisjointPartitioningDataHandle.
  /**
   * Creates a new DisjointPartitioningDataHandle and by default initializes the
   * result vector with the current MPI rank. If you have already done that
   * externally, you can skip the initialization.
   *
   * \param space_         The GridFunctionSpace to operate on.
   * \param v_           The result vector.
   * \param init_vector  Flag to control whether the result vector will be initialized.
   */
  DisjointPartitioningDataHandle(const SpaceType& space_, VectorType& v_, bool init_vector = true)
    : BaseType(space_,
               v_,
               DisjointPartitioningGatherScatter<typename VectorType::ScalarType>(space_.grid_layer().comm().rank()))
  {
    if (init_vector)
      v_.set_all(space_.grid_layer().comm().rank());
  }
};


//! GatherScatter functor for marking shared DOFs.
/**
 * This functor will mark all DOFs that exist on multiple processes.
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::All_All_Interface and the result vector must be initialized with false.
 */
struct SharedDOFGatherScatter
{

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool gather(MessageBuffer& buff, const Entity& e, LocalView& local_view) const
  {
    buff.write(local_view.size() > 0);
    return false;
  }

  template <typename MessageBuffer, typename Entity, typename LocalView>
  bool scatter(MessageBuffer& buff, std::size_t n, const Entity& e, LocalView& local_view) const
  {
    bool remote_entity_has_dofs;
    buff.read(remote_entity_has_dofs);

    for (std::size_t i = 0; i < local_view.size(); ++i) {
      local_view[i] |= remote_entity_has_dofs;
    }
    return true;
  }
};


//! Data handle for marking shared DOFs.
/**
 * This data handle will mark all DOFs that exist on multiple processes.
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::All_All_Interface and the result vector must be initialized with false.
 */
template <class SpaceType, class VectorType>
class SharedDOFDataHandle
    : public SpaceDataHandle<SpaceType, VectorType, SharedDOFGatherScatter, EntityDataCommunicationDescriptor<bool>>
{
  typedef SpaceDataHandle<SpaceType, VectorType, SharedDOFGatherScatter, EntityDataCommunicationDescriptor<bool>>
      BaseType;

  //  static_assert((std::is_same<typename VectorType::ScalarType, bool>::value),
  //                "SharedDOFDataHandle expects a vector of bool values");

public:
  //! Creates a new SharedDOFDataHandle.
  /**
   * Creates a new SharedDOFDataHandle and by default initializes the result vector
   * with the correct value of false. If you have already done that externally,
   * you can skip the initialization.
   *
   * \param space_         The GridFunctionSpace to operate on.
   * \param v_           The result vector.
   * \param init_vector  Flag to control whether the result vector will be initialized.
   */
  SharedDOFDataHandle(const SpaceType& space_, VectorType& v_, bool init_vector = true)
    : BaseType(space_, v_)
  {
    if (init_vector)
      v_.set_all(false);
  }
};


//! Data handle for collecting set of neighboring MPI ranks.
/**
 * This data handle collects the MPI ranks of all processes that share grid entities
 * with attached DOFs.
 *
 * \note In order to work correctly, the data handle must be communicated on the
 * Dune::All_All_Interface.
 */
template <typename SpaceType, typename RankIndex>
class SpaceNeighborDataHandle : public Dune::CommDataHandleIF<SpaceNeighborDataHandle<SpaceType, RankIndex>, RankIndex>
{

  // We deliberately avoid using the SpaceDataHandle here, as we don't want to incur the
  // overhead of invoking the whole SpaceType infrastructure.

public:
  typedef RankIndex DataType;

  SpaceNeighborDataHandle(const SpaceType& space, RankIndex rank, std::set<RankIndex>& neighbors)
    : space_(space)
    , _rank(rank)
    , _neighbors(neighbors)
  {
  }

  bool contains(int dim, int codim) const
  {
    return dim == SpaceType::dimDomain && codim == 0;
  }

  bool fixedsize(int /*dim*/, int /*codim*/) const
  {
    // We always send a single value, the MPI rank.
    return true;
  }

  template <typename Entity>
  size_t size(Entity& /*e*/) const
  {
    return 1;
  }

  template <typename MessageBuffer, typename Entity>
  void gather(MessageBuffer& buff, const Entity& /*e*/) const
  {
    buff.write(_rank);
  }

  template <typename MessageBuffer, typename Entity>
  void scatter(MessageBuffer& buff, const Entity& e, size_t n)
  {
    RankIndex rank;
    buff.read(rank);
    _neighbors.insert(rank);
  }

private:
  const SpaceType& space_;
  const RankIndex _rank;
  std::set<RankIndex>& _neighbors;
};


} // namespace GDT
} // namespace Dune
#endif // DUNE_GDT_DATAHANDLES_HH
