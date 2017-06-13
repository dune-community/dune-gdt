#ifndef DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH
#define DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH

#include <dune/grid/common/datahandleif.hh>

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType>
class DiscreteFunctionDataHandle
    : public Dune::CommDataHandleIF<DiscreteFunctionDataHandle<DiscreteFunctionType>,
                                    typename DiscreteFunctionType::SpaceType::RangeFieldType>
{
public:
  DiscreteFunctionDataHandle(DiscreteFunctionType& discrete_function, bool fixed_size = true)
    : discrete_function_(discrete_function)
    , fixed_size_(fixed_size)
  {
  }

  //! export type of data for message buffer
  typedef typename DiscreteFunctionType::SpaceType::DomainFieldType DataType;

  //! returns true if data for this codim should be communicated
  bool contains(int /*dim*/, int codim) const
  {
    return (codim == 0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize(int /*dim*/, int /*codim*/) const
  {
    return fixed_size_;
  }

  /*! how many objects of type DataType have to be sent for a given entity

     Note: Only the sender side needs to know this size.
   */
  template <class EntityType>
  std::enable_if_t<EntityType::codimension != 0, size_t> size(const EntityType& /*entity*/) const
  {
    return 0;
  }

  /*! how many objects of type DataType have to be sent for a given entity
     Note: Only the sender side needs to know this size.
   */
  template <class EntityType>
  std::enable_if_t<EntityType::codimension == 0, size_t> size(const EntityType& entity) const
  {
    return discrete_function_.space().mapper().numDofs(entity);
  }

  //! pack data from user to message buffer
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0> gather(MessageBuffer& /*buff*/, const EntityType& /*entity*/) const
  {
  }

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0> gather(MessageBuffer& buff, const EntityType& entity) const
  {
    const auto& mapper = discrete_function_.space().mapper();
    const auto& vector = discrete_function_.vector();
    for (size_t ii = 0; ii < mapper.numDofs(entity); ++ii)
      buff.write(vector.get_entry(mapper.mapToGlobal(entity, ii)));
  }

  /*! unpack data from message buffer to user
     n is the number of objects sent by the sender
   */
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0>
  scatter(MessageBuffer& /*buff*/, const EntityType& /*entity*/, size_t /*n*/)
  {
  }

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0> scatter(MessageBuffer& buff, const EntityType& entity, size_t n)
  {
    const auto& mapper = discrete_function_.space().mapper();
    auto& vector = discrete_function_.vector();
    assert(mapper.numDofs(entity) == n);
    // we need this intermediate double because we cannot get a reference to an entry in vector
    typename DiscreteFunctionType::SpaceType::RangeFieldType entry;
    for (size_t ii = 0; ii < n; ++ii) {
      buff.read(entry);
      vector.set_entry(mapper.mapToGlobal(entity, ii), entry);
    }
  }

private:
  DiscreteFunctionType& discrete_function_;
  const bool fixed_size_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH
