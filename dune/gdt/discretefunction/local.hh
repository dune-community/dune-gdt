// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH

#include <vector>
#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/mapper/interface.hh>

namespace Dune {
namespace GDT {


template< class VectorImp >
class ConstLocalDoFVector
{
  static_assert(std::is_base_of
                < Stuff::LA::VectorInterface< typename VectorImp::Traits, typename VectorImp::Traits::ScalarType >,
                  VectorImp >::value,
                "VectorImp has to be derived from Stuff::LA::VectorInterface!");
public:
  typedef VectorImp VectorType;
  typedef typename VectorType::ScalarType ScalarType;

  template< class M, class EntityType >
  ConstLocalDoFVector(const MapperInterface< M >& mapper,
                      const EntityType& entity,
                      const VectorType& vector)
    : vector_(vector)
    , indices_(mapper.numDofs(entity))
  {
    mapper.globalIndices(entity, indices_);
  }

  ~ConstLocalDoFVector() {}

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
  Dune::DynamicVector< size_t > indices_;

private:
  template< class V >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ConstLocalDoFVector< V >& /*vector*/);
}; // class ConstLocalDoFVector


template< class V >
std::ostream& operator<<(std::ostream& out, const ConstLocalDoFVector< V >& vector)
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


template< class VectorImp >
class LocalDoFVector
  : public ConstLocalDoFVector< VectorImp >
{
  typedef ConstLocalDoFVector< VectorImp > BaseType;
public:
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::ScalarType ScalarType;

  template< class M, class EntityType >
  LocalDoFVector(const MapperInterface< M >& mapper,
                 const EntityType& entity,
                 VectorType& vector)
    : BaseType(mapper, entity, vector)
    , vector_(vector)
  {}

  ~LocalDoFVector() {}

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

private:
  using BaseType::indices_;
  VectorType& vector_;
}; // class LocalDoFVector


template< class SpaceImp, class VectorImp >
class ConstLocalDiscreteFunction
  : public Stuff::LocalfunctionInterface< typename SpaceImp::EntityType,
                                          typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
                                          typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
{
  static_assert(std::is_base_of< SpaceInterface< typename SpaceImp::Traits,
                                                 SpaceImp::dimDomain,
                                                 SpaceImp::dimRange,
                                                 SpaceImp::dimRangeCols >,
                                 SpaceImp >::value,
                "SpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of
                < Dune::Stuff::LA::VectorInterface< typename VectorImp::Traits, typename VectorImp::Traits::ScalarType >,
                  VectorImp >::value,
                "VectorImp has to be derived from Stuff::LA::VectorInterface!");
  static_assert(std::is_same< typename SpaceImp::RangeFieldType, typename VectorImp::ScalarType >::value,
                "Types do not match!");
  typedef Stuff::LocalfunctionInterface
      < typename SpaceImp::EntityType, typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
        typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
    BaseType;
  typedef ConstLocalDiscreteFunction< SpaceImp, VectorImp > ThisType;
public:
  typedef SpaceImp                      SpaceType;
  typedef VectorImp                     VectorType;
  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType       DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRangeRows = BaseType::dimRangeCols;
  static const unsigned int                 dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType      RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;
private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;

public:

  ConstLocalDiscreteFunction(const SpaceType& space, const VectorType& global_vector, const EntityType& ent)
    : BaseType(ent)
    , space_(space)
    , base_(new BaseFunctionSetType(space_.base_function_set(this->entity())))
    , global_indices_(space_.mapper().globalIndices(entity_))
    , const_local_DoF_vector_(map_to_local(global_vector, global_indices_))
  {
    assert(const_local_DoF_vector_.size() == base_->size());
  }

  ConstLocalDiscreteFunction(ThisType&& source) = default;

  ConstLocalDiscreteFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~ConstLocalDiscreteFunction() {}

  const BaseFunctionSetType& base() const
  {
    return *base_;
  }

  const VectorType& vector() const
  {
    return const_local_DoF_vector_;
  }

  virtual size_t order() const override
  {
    return base_->order();
  }

  virtual void evaluate(const DomainType& xx, RangeType& ret) const override
  {
    assert(this->is_a_valid_point(xx));
    ret *= 0.0;
    std::vector<RangeType> tmpBaseValues(base_->size(), RangeType(0));
    assert(localVector_->size() == tmpBaseValues.size());
    base_->evaluate(xx, tmpBaseValues);
    for (size_t ii = 0; ii < const_local_DoF_vector_.size(); ++ii) {
      tmpBaseValues[ii] *= const_local_DoF_vector_[ii];
      ret += tmpBaseValues[ii];
    }
  } // ... evaluate(...)

  virtual void jacobian(const DomainType& xx, JacobianRangeType& ret) const override
  {
    assert(this->is_a_valid_point(xx));
    ret *= RangeFieldType(0);
    std::vector<JacobianRangeType> tmpBaseJacobianValues(base_->size(), JacobianRangeType(0));
    assert(localVector_->size() == tmpBaseJacobianValues.size());
    base_->jacobian(xx, tmpBaseJacobianValues);
    for (size_t ii = 0; ii < const_local_DoF_vector_.size(); ++ii) {
      tmpBaseJacobianValues[ii] *= const_local_DoF_vector_[ii];
      ret += tmpBaseJacobianValues[ii];
    }
  } // ... jacobian(...)

  using BaseType::evaluate;
  using BaseType::jacobian;

protected:
  static VectorType map_to_local(const VectorType& global_vector, const DynamicVector< size_t >& indices)
  {
    VectorType ret(indices.size());
    for (size_t ii = 0; ii < indices.size(); ++ii) {
      assert(indices[ii] < global_vector.size());
      ret[ii] = global_vector[indices[ii]];
    }
    return ret;
  } // ... map_to_local(...)

  using BaseType::entity_;
  const SpaceType& space_;
  const std::unique_ptr< const BaseFunctionSetType > base_;
  const DynamicVector< size_t > global_indices_;
  const VectorType const_local_DoF_vector_;
}; // class ConstLocalDiscreteFunction


template< class SpaceImp, class VectorImp >
class LocalDiscreteFunction
  : public ConstLocalDiscreteFunction< SpaceImp, VectorImp >
{
  typedef ConstLocalDiscreteFunction< SpaceImp, VectorImp > BaseType;
  typedef LocalDiscreteFunction< SpaceImp, VectorImp >      ThisType;
public:
  typedef typename BaseType::SpaceType  SpaceType;
  typedef typename BaseType::VectorType VectorType;
  typedef typename BaseType::EntityType EntityType;

  typedef typename BaseType::DomainFieldType  DomainFieldType;
  static const unsigned int                   dimDomain = BaseType::dimDomain;
  typedef typename BaseType::DomainType       DomainType;

  typedef typename BaseType::RangeFieldType RangeFieldType;
  static const unsigned int                 dimRangeRows = BaseType::dimRangeCols;
  static const unsigned int                 dimRangeCols = BaseType::dimRangeCols;
  typedef typename BaseType::RangeType      RangeType;

  typedef typename BaseType::JacobianRangeType JacobianRangeType;
private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;

public:
  LocalDiscreteFunction(const SpaceType& space, VectorType& global_vector, const EntityType& ent)
    : BaseType(space, global_vector, ent)
    , global_vector_(global_vector)
    , local_DoF_vector_(const_local_DoF_vector_)
  {}

  //! previous comment questioned validity, defaulting this doesn't touch that question
  LocalDiscreteFunction(ThisType&& source) = default;

  LocalDiscreteFunction(const ThisType& other) = delete;

  ThisType& operator=(const ThisType& other) = delete;

  virtual ~LocalDiscreteFunction()
  {
    for (size_t ii = 0; ii < global_indices_.size(); ++ii) {
      assert(global_indices_[ii] < global_vector_.size());
      global_vector_[global_indices_[ii]] = local_DoF_vector_[ii];
    }
  } // ~LocalDiscreteFunction(...)

  VectorType& vector()
  {
    return local_DoF_vector_;
  }

private:
  using BaseType::space_;
  using BaseType::entity_;
  using BaseType::base_;
  using BaseType::global_indices_;
  using BaseType::const_local_DoF_vector_;
  VectorType& global_vector_;
  VectorType local_DoF_vector_;
}; // class LocalDiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
