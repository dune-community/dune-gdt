// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_DATAHANDLE_HH
#define DUNE_GDT_OPERATORS_FV_DATAHANDLE_HH

#include <dune/grid/common/datahandleif.hh>

namespace Dune {
namespace GDT {


template <class ReconstructedValuesType, class GridLayerType>
class ReconstructionDataHandle
    : public Dune::CommDataHandleIF<ReconstructionDataHandle<ReconstructedValuesType, GridLayerType>,
                                    typename ReconstructedValuesType::value_type::mapped_type::value_type>
{
  typedef typename ReconstructedValuesType::value_type::mapped_type RangeType;

public:
  ReconstructionDataHandle(ReconstructedValuesType& reconstructed_values, const GridLayerType& grid_layer)
    : reconstructed_values_(reconstructed_values)
    , grid_layer_(grid_layer)
  {
  }

  //! export type of data for message buffer
  typedef typename RangeType::value_type DataType;

  //! returns true if data for this codim should be communicated
  bool contains(int /*dim*/, int codim) const
  {
    return (codim == 0);
  }

  //! returns true if size per entity of given dim and codim is a constant
  bool fixedsize(int /*dim*/, int /*codim*/) const
  {
    return true;
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
    return reconstructed_values_[grid_layer_.indexSet().index(entity)].size() * RangeType::dimension;
  }

  //! pack data from user to message buffer
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0> gather(MessageBuffer& /*buff*/, const EntityType& /*entity*/) const
  {
  }

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0> gather(MessageBuffer& buff, const EntityType& entity) const
  {
    for (const auto& pair : reconstructed_values_[grid_layer_.indexSet().index(entity)])
      for (size_t rr = 0; rr < RangeType::dimension; ++rr)
        buff.write(pair.second[rr]);
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
  std::enable_if_t<EntityType::codimension == 0> scatter(MessageBuffer& buff, const EntityType& entity, size_t /*n*/)
  {
    for (auto& pair : reconstructed_values_[grid_layer_.indexSet().index(entity)])
      for (size_t rr = 0; rr < RangeType::dimension; ++rr)
        buff.read(pair.second[rr]);
  }

private:
  ReconstructedValuesType& reconstructed_values_;
  const GridLayerType& grid_layer_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_DATAHANDLE_HH
