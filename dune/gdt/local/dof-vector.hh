// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2018)
//   René Fritze     (2016, 2018)
//   René Milk       (2017)

#ifndef DUNE_GDT_LOCAL_DOF_VECTOR_HH
#define DUNE_GDT_LOCAL_DOF_VECTOR_HH

#include <dune/xt/common/vector.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/bound-object.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {
namespace internal {


// forward, required for the default in the declaration of ConstLocalDofVector
template <class Vector, class GridView>
class ConstLocalDofVectorTraits;


} // namespace internal


// forwards, required for the traits
template <class Vector, class GridView, class Traits = internal::ConstLocalDofVectorTraits<Vector, GridView>>
class ConstLocalDofVector;

template <class Vector, class GridView>
class LocalDofVector;

// required to avoid cyclic inclusion
template <class GV>
class MapperInterface;


namespace internal {


template <class Vector, class GridView>
class ConstLocalDofVectorTraits
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

public:
  using derived_type = ConstLocalDofVector<Vector, GridView>;
  using ScalarType = typename Vector::ScalarType;
  using RealType = typename Vector::RealType;
};


template <class Vector, class GridView>
class LocalDofVectorTraits
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

public:
  using derived_type = LocalDofVector<Vector, GridView>;
  using ScalarType = typename Vector::ScalarType;
  using RealType = typename Vector::RealType;
};


} // namespace internal


/**
 * \note Since all implementations of XT::LA::VectorInterface are assumed to be thread safe, no special care needs to be
 *       taken here.
 */
template <class Vector, class GridView, class Traits>
class ConstLocalDofVector
  : public XT::LA::VectorInterface<Traits>
  , public XT::Grid::ElementBoundObject<XT::Grid::extract_entity_t<GridView>>
{

  using ThisType = ConstLocalDofVector<Vector, GridView, Traits>;
  using BaseType = XT::LA::VectorInterface<Traits>;

public:
  using typename BaseType::ScalarType;
  using VectorType = Vector;
  using MapperType = MapperInterface<GridView>;
  using typename XT::Grid::ElementBoundObject<XT::Grid::extract_entity_t<GridView>>::ElementType;

  ConstLocalDofVector(const MapperType& mapper, const VectorType& global_vector)
    : mapper_(mapper)
    , global_vector_(global_vector)
    , global_DoF_indices_(mapper_.max_local_size())
    , size_(0)
  {}

  ConstLocalDofVector(const ThisType& other) = default;
  ConstLocalDofVector(ThisType&& source) = default;

protected:
  void post_bind(const ElementType& ele) override final
  {
    mapper_.global_indices(ele, global_DoF_indices_);
    size_ = mapper_.local_size(ele);
    DUNE_THROW_IF(global_DoF_indices_.size() < size_,
                  Exceptions::dof_vector_error,
                  "This must not happen, the mapper is broken!");
  }

public:
  size_t size() const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return size_;
  }

  void resize(const size_t new_size)
  {
    DUNE_THROW_IF(this->size() != new_size, Exceptions::dof_vector_error, "this does not make sense!");
  }

  void add_to_entry(const size_t /*ii*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(Exceptions::dof_vector_error, "a ConstLocalDofVector is not mutable!");
  }

  void set_entry(const size_t /*ii*/, const ScalarType& /*value*/)
  {
    DUNE_THROW(Exceptions::dof_vector_error, "a ConstLocalDofVector is not mutable!");
  }

  ScalarType get_entry(const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_.get_entry(global_DoF_indices_[ii]);
  }

protected:
  ScalarType& get_unchecked_ref(const size_t /*ii*/)
  {
    DUNE_THROW(Exceptions::dof_vector_error, "a ConstLocalDofVector is not mutable!");
    return dummy_scalar_to_silence_the_warning_;
  }

  const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    // This is not optimal, but global_vector_.get_unchecked_ref is protected.
    return global_vector_[global_DoF_indices_[ii]];
  }

public:
  ScalarType& operator[](const size_t /*ii*/)
  {
    DUNE_THROW(Exceptions::dof_vector_error, "a ConstLocalDofVector is not mutable!");
    return dummy_scalar_to_silence_the_warning_;
  }

  const ScalarType& operator[](const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_[global_DoF_indices_[ii]];
  }

protected:
  friend BaseType; // To allow access to get_unchecked_ref().

  const MapperType& mapper_;

private:
  const VectorType& global_vector_;

protected:
  DynamicVector<size_t> global_DoF_indices_;
  size_t size_;

private:
  ScalarType dummy_scalar_to_silence_the_warning_;
}; // class ConstLocalDofVector


template <class Vector, class GridView>
class LocalDofVector : public ConstLocalDofVector<Vector, GridView, internal::LocalDofVectorTraits<Vector, GridView>>
{

  using ThisType = LocalDofVector<Vector, GridView>;
  using BaseType = ConstLocalDofVector<Vector, GridView, internal::LocalDofVectorTraits<Vector, GridView>>;

public:
  using typename BaseType::ScalarType;
  using VectorType = Vector;
  using MapperType = MapperInterface<GridView>;

  LocalDofVector(const MapperType& mapper, VectorType& global_vector)
    : BaseType(mapper, global_vector)
    , global_vector_(global_vector)
  {}

  LocalDofVector(const ThisType&) = default;
  LocalDofVector(ThisType&&) = default;

  void add_to_entry(const size_t ii, const ScalarType& value)
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    global_vector_.add_to_entry(global_DoF_indices_[ii], value);
  }

  void set_entry(const size_t ii, const ScalarType& value)
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    global_vector_.set_entry(global_DoF_indices_[ii], value);
  }

  ScalarType get_entry(const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_.get_entry(global_DoF_indices_[ii]);
  }

protected:
  ScalarType& get_unchecked_ref(const size_t ii)
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_[global_DoF_indices_[ii]];
  }

  const ScalarType& get_unchecked_ref(const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_[global_DoF_indices_[ii]];
  }

public:
  ScalarType& operator[](const size_t ii)
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_[global_DoF_indices_[ii]];
  }

  const ScalarType& operator[](const size_t ii) const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    assert(ii < size_);
    return global_vector_[global_DoF_indices_[ii]];
  }

private:
  // To allow access to get_unchecked_ref().
  friend XT::LA::VectorInterface<internal::LocalDofVectorTraits<Vector, GridView>>;

  using BaseType::mapper_;
  VectorType& global_vector_;
  using BaseType::global_DoF_indices_;
  using BaseType::size_;
}; // class LocalDofVector


} // namespace GDT
namespace XT {
namespace Common {


template <class Vector, class GridView, class Traits>
struct VectorAbstraction<GDT::ConstLocalDofVector<Vector, GridView, Traits>>
  : public LA::internal::VectorAbstractionBase<GDT::ConstLocalDofVector<Vector, GridView, Traits>>
{};


template <class Vector, class GridView>
struct VectorAbstraction<GDT::LocalDofVector<Vector, GridView>>
  : public LA::internal::VectorAbstractionBase<GDT::LocalDofVector<Vector, GridView>>
{};


} // namespace Common
} // namespace XT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_DOF_VECTOR_HH
