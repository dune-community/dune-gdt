// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH
#define DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH

#include <dune/grid/common/datahandleif.hh>

#include <dune/xt/common/unused.hh>

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType>
class DiscreteFunctionDataHandle
  : public Dune::CommDataHandleIF<DiscreteFunctionDataHandle<DiscreteFunctionType>,
                                  typename DiscreteFunctionType::SpaceType::R>
{
public:
  DiscreteFunctionDataHandle(DiscreteFunctionType& discrete_function, bool fixed_size = true)
    : discrete_function_(discrete_function)
    , mapper_(discrete_function_.space().mapper())
    , vector_(discrete_function_.dofs().vector())
    , fixed_size_(fixed_size)
  {}

  //! export type of data for message buffer
  typedef typename DiscreteFunctionType::SpaceType::D DataType;

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
    return discrete_function_.space().mapper().local_size(entity);
  }

  //! pack data from user to message buffer
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0> gather(MessageBuffer& /*buff*/, const EntityType& /*entity*/) const
  {}

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0> gather(MessageBuffer& buff, const EntityType& entity) const
  {
    const auto global_indices = mapper_.global_indices(entity);
    for (const auto& index : global_indices)
      buff.write(vector_.get_entry(index));
  }

  /*! unpack data from message buffer to user
     n is the number of objects sent by the sender
   */
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0>
  scatter(MessageBuffer& /*buff*/, const EntityType& /*entity*/, size_t /*n*/)
  {}

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0>
  scatter(MessageBuffer& buff, const EntityType& entity, size_t DXTC_DEBUG_ONLY(n))
  {
    assert(mapper_.local_size(entity) == n);
    const auto global_indices = mapper_.global_indices(entity);
    for (const auto& index : global_indices)
      buff.read(vector_[index]);
  }

private:
  DiscreteFunctionType& discrete_function_;
  const typename DiscreteFunctionType::SpaceType::MapperType& mapper_;
  typename DiscreteFunctionType::VectorType& vector_;
  const bool fixed_size_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DATAHANDLE_HH
