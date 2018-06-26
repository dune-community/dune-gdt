// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_ASSEMBLER_CODIM0MATRIXDATAHANDLE_HH
#define DUNE_GDT_ASSEMBLER_CODIM0MATRIXDATAHANDLE_HH

#include <cmath>

#include <dune/grid/common/datahandleif.hh>

#include <dune/xt/common/unused.hh>

namespace Dune {
namespace GDT {


// datahandle for a Dune::XT::LA::MatrixInterface matrix that can be assembled from local matrices on each entity
template <class MatrixType, class SpaceType>
class MatrixDataHandle
    : public Dune::CommDataHandleIF<MatrixDataHandle<MatrixType, SpaceType>, typename MatrixType::ScalarType>
{
public:
  MatrixDataHandle(MatrixType& matrix, const SpaceType& space, bool fixed_size = true)
    : matrix_(matrix)
    , space_(space)
    , fixed_size_(fixed_size)
  {
  }

  //! export type of data for message buffer
  typedef typename MatrixType::ScalarType DataType;

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
    return std::pow(space_.mapper().numDofs(entity), 2);
  }

  //! pack data from user to message buffer
  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension != 0> gather(MessageBuffer& /*buff*/, const EntityType& /*entity*/) const
  {
  }

  template <class MessageBuffer, class EntityType>
  std::enable_if_t<EntityType::codimension == 0> gather(MessageBuffer& buff, const EntityType& entity) const
  {
    const auto& mapper = space_.mapper();
    const auto& num_local_dofs = mapper.numDofs(entity);
    for (size_t ii = 0; ii < num_local_dofs; ++ii)
      for (size_t jj = 0; jj < num_local_dofs; ++jj)
        buff.write(matrix_.get_entry(mapper.mapToGlobal(entity, ii), mapper.mapToGlobal(entity, jj)));
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
  std::enable_if_t<EntityType::codimension == 0>
  scatter(MessageBuffer& buff, const EntityType& entity, size_t DXTC_DEBUG_ONLY(n))
  {
    const auto& mapper = space_.mapper();
    const auto& num_local_dofs = mapper.numDofs(entity);
    assert(num_local_dofs * num_local_dofs == n);
    for (size_t ii = 0; ii < num_local_dofs; ++ii)
      for (size_t jj = 0; jj < num_local_dofs; ++jj)
        buff.read(matrix_.get_entry_ref(mapper.mapToGlobal(entity, ii), mapper.mapToGlobal(entity, jj)));
  }

private:
  MatrixType& matrix_;
  const SpaceType& space_;
  const bool fixed_size_;
};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_CODIM0MATRIXDATAHANDLE_HH
