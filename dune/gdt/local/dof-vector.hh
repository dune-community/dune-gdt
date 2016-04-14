// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)

#ifndef DUNE_GDT_LOCAL_DOF_VECTOR_HH
#define DUNE_GDT_LOCAL_DOF_VECTOR_HH

#include <ostream>

#include <dune/stuff/la/container/vector-interface.hh>

#include <dune/gdt/spaces/mapper/interfaces.hh>

namespace Dune {
namespace GDT {


template <class VectorImp>
class ConstLocalDoFVector
{
  static_assert(
      std::is_base_of<Stuff::LA::VectorInterface<typename VectorImp::Traits, typename VectorImp::Traits::ScalarType>,
                      VectorImp>::value,
      "VectorImp has to be derived from Stuff::LA::VectorInterface!");

public:
  typedef VectorImp VectorType;
  typedef typename VectorType::ScalarType ScalarType;

  template <class M, class EntityType>
  ConstLocalDoFVector(const MapperInterface<M>& mapper, const EntityType& entity, const VectorType& vector)
    : vector_(vector)
    , indices_(mapper.numDofs(entity))
  {
    mapper.globalIndices(entity, indices_);
  }

  ~ConstLocalDoFVector()
  {
  }

  size_t size() const
  {
    return indices_.size();
  }

  ScalarType get(const size_t ii) const
  {
    assert(ii < indices_.size());
    return vector_.get_entry(indices_[ii]);
  }

private:
  const VectorType& vector_;

protected:
  Dune::DynamicVector<size_t> indices_;

private:
  template <class V>
  friend std::ostream& operator<<(std::ostream& /*out*/, const ConstLocalDoFVector<V>& /*vector*/);
}; // class ConstLocalDoFVector


template <class V>
std::ostream& operator<<(std::ostream& out, const ConstLocalDoFVector<V>& vector)
{
  out << "[";
  const size_t sz = vector.size();
  if (sz > 0) {
    out << vector.get(0);
    for (size_t ii = 1; ii < sz; ++ii)
      out << ", " << vector.get(ii);
  } else
    out << " ";
  out << "]";
  return out;
} // ... operator<<(...)


template <class VectorImp>
class LocalDoFVector : public ConstLocalDoFVector<VectorImp>
{
  typedef ConstLocalDoFVector<VectorImp> BaseType;

public:
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ScalarType ScalarType;

  template <class M, class EntityType>
  LocalDoFVector(const MapperInterface<M>& mapper, const EntityType& entity, VectorType& vector)
    : BaseType(mapper, entity, vector)
    , vector_(vector)
  {
  }

  ~LocalDoFVector()
  {
  }

  void set(const size_t ii, const ScalarType& val)
  {
    assert(ii < indices_.size());
    vector_.set_entry(indices_[ii], val);
  }

  void add(const size_t ii, const ScalarType& val)
  {
    assert(ii < indices_.size());
    vector_.add_to_entry(indices_[ii], val);
  }

  template <class OtherVectorImp>
  void add(const OtherVectorImp& vector)
  {
    assert(vector.size() == indices_.size());
    for (size_t ii = 0; ii < indices_.size(); ++ii)
      add(ii, vector[ii]);
  }

private:
  using BaseType::indices_;
  VectorType& vector_;
}; // class LocalDoFVector


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_DOF_VECTOR_HH
