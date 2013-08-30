// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
#define DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH

#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interface.hh>
#include <dune/stuff/common/vector.hh>
#include <dune/stuff/common/matrix.hh>

#include <dune/gdt/mapper/interface.hh>

namespace Dune {
namespace GDT {


template< class VectorImp >
class LocalDoFVector
{
public:
  typedef typename Dune::Stuff::LA::VectorInterface< typename VectorImp::Traits >::derived_type VectorType;
  typedef typename VectorType::ElementType ElementType;

  template< class M, class EntityType >
  LocalDoFVector(const MapperInterface< M >& mapper,
                 const EntityType& entity,
                 VectorType& vector)
    : indices_(mapper.numDofs(entity))
    , vector_(vector)
  {
    mapper.globalIndices(entity, indices_);
  }

  size_t size() const
  {
    return indices_.size();
  }

  ElementType get(const size_t ii) const
  {
    assert(ii < indices_.size());
    return vector_.get(indices_[ii]);
  }

  void set(const size_t ii, const ElementType& val)
  {
    assert(ii < indices_.size());
    vector_.set(indices_[ii], val);
  }

  void add(const size_t ii, const ElementType& val)
  {
    assert(ii < indices_.size());
    vector_.add(indices_[ii], val);
  }

private:
  Dune::DynamicVector< size_t > indices_;
  VectorType& vector_;
}; // class LocalDoFVector


// forward, includes are below
template< class SpaceImp, class VectorImp >
class DiscreteFunctionDefaultConst;
template< class SpaceImp, class VectorImp >
class DiscreteFunctionDefault;


// forward, to be used in the traits
template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocalConst;


template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocalConstTraits
{
public:
  typedef DiscreteFunctionLocalConst< SpaceImp, VectorImp > derived_type;
  typedef typename DiscreteFunctionDefaultConst< SpaceImp, VectorImp >::EntityType EntityType;
};


template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocalConst
  : public Dune::Stuff::LocalFunctionInterface< DiscreteFunctionLocalConstTraits< SpaceImp, VectorImp >,
                                                typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
                                                typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
{
public:
  typedef DiscreteFunctionLocalConstTraits< SpaceImp, VectorImp > Traits;
  typedef DiscreteFunctionDefaultConst< SpaceImp, VectorImp >         DiscreteFunctionType;
  typedef typename DiscreteFunctionType::SpaceType                    SpaceType;
  typedef LocalDoFVector< typename DiscreteFunctionType::VectorType > LocalDoFVectorType;
  typedef typename DiscreteFunctionType::EntityType                   EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int                   dimRangeCols = SpaceType::dimRangeCols;
private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
public:
  typedef typename BaseFunctionSetType::DomainType        DomainType;
  typedef typename BaseFunctionSetType::RangeType         RangeType;
  typedef typename BaseFunctionSetType::JacobianRangeType JacobianRangeType;

  DiscreteFunctionLocalConst(const DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : function_(discreteFunction)
    , entity_(entity)
    , base_(function_.space().baseFunctionSet(entity_))
    , localVector_(function_.space().mapper(), entity_, const_cast< typename DiscreteFunctionType::VectorType& >(*(function_.vector())))
    , tmpBaseValues_(base_.size(), RangeType(0))
    , tmpBaseJacobianValues_(base_.size(), JacobianRangeType(0))
  {
    assert(localVector_.size() == base_.size());
  }

  const DiscreteFunctionType& discreteFunction() const
  {
    return function_;
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const LocalDoFVectorType& vector() const
  {
    return localVector_;
  }

  int order() const
  {
    return base_.order();
  }

  size_t size() const
  {
    return base_.size();
  }

  void evaluate(const DomainType& x, RangeType& ret) const
  {
    Dune::Stuff::Common::clear(ret);
    assert(localVector_.size() == tmpBaseValues_.size());
    base_.evaluate(x, tmpBaseValues_);
    for (size_t ii = 0; ii < localVector_.size(); ++ii) {
      tmpBaseValues_[ii] *= localVector_.get(ii);
      ret += tmpBaseValues_[ii];
    }
  } // ... evaluate(...)

  RangeType evaluate(const DomainType& xx) const
  {
    RangeType ret;
    evaluate(xx, ret);
    return ret;
  }

  void jacobian(const DomainType& x, JacobianRangeType& ret) const
  {
    Dune::Stuff::Common::clear(ret);
    assert(localVector_.size() == tmpBaseJacobianValues_.size());
    base_.jacobian(x, tmpBaseJacobianValues_);
    for (size_t ii = 0; ii < localVector_.size(); ++ii) {
      tmpBaseJacobianValues_[ii] *= localVector_.get(ii);
      ret += tmpBaseJacobianValues_[ii];
    }
  } // ... jacobian(...)

  JacobianRangeType jacobian(const DomainType& xx) const
  {
    JacobianRangeType ret;
    jacobian(xx, ret);
    return ret;
  }

private:
  const DiscreteFunctionType& function_;
  const EntityType& entity_;
  const BaseFunctionSetType base_;
protected:
  LocalDoFVectorType localVector_;
private:
  mutable std::vector< RangeType > tmpBaseValues_;
  mutable std::vector< JacobianRangeType > tmpBaseJacobianValues_;
}; // class DiscreteFunctionLocalConst


// forward, to be used in the traits
template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocal;


template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocalTraits
{
public:
  typedef DiscreteFunctionLocal< SpaceImp, VectorImp > derived_type;
  typedef typename DiscreteFunctionDefault< SpaceImp, VectorImp >::EntityType EntityType;
};


template< class SpaceImp, class VectorImp >
class DiscreteFunctionLocal
  : public DiscreteFunctionLocalConst< SpaceImp, VectorImp >
  , public Dune::Stuff::LocalFunctionInterface< DiscreteFunctionLocalTraits< SpaceImp, VectorImp >,
                                                typename SpaceImp::DomainFieldType, SpaceImp::dimDomain,
                                                typename SpaceImp::RangeFieldType, SpaceImp::dimRange, SpaceImp::dimRangeCols >
{
  typedef DiscreteFunctionLocalConst< SpaceImp, VectorImp > BaseType;
public:
  typedef DiscreteFunctionDefault< SpaceImp, VectorImp >              DiscreteFunctionType;
  typedef typename DiscreteFunctionType::SpaceType                    SpaceType;
  typedef LocalDoFVector< typename DiscreteFunctionType::VectorType > LocalDoFVectorType;
  typedef typename DiscreteFunctionType::EntityType                   EntityType;
  typedef typename SpaceType::DomainFieldType DomainFieldType;
  static const unsigned int                   dimDomain = SpaceType::dimDomain;
  typedef typename SpaceType::RangeFieldType  RangeFieldType;
  static const unsigned int                   dimRange = SpaceType::dimRange;
  static const unsigned int                   dimRangeCols = SpaceType::dimRangeCols;
private:
  typedef typename SpaceType::BaseFunctionSetType BaseFunctionSetType;
public:
  typedef typename BaseFunctionSetType::DomainType        DomainType;
  typedef typename BaseFunctionSetType::RangeType         RangeType;
  typedef typename BaseFunctionSetType::JacobianRangeType JacobianRangeType;

  DiscreteFunctionLocal(const DiscreteFunctionType& discreteFunction, const EntityType& entity)
    : BaseType(discreteFunction, entity)
  {}

  LocalDoFVectorType& vector()
  {
    return BaseType::localVector_;
  }
}; // class DiscreteFunctionLocal


} // namespace GDT
} // namespace Dune

#include "default.hh"

#endif // DUNE_GDT_DISCRETEFUNCTION_LOCAL_HH
