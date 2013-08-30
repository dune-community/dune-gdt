// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_DISCRETEFUNCTION_DIFFERENCE_HH
#define DUNE_GDT_DISCRETEFUNCTION_DIFFERENCE_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {

namespace LocalFunction {

template< class EntityType, class MinuendType, class SubtrahendType  >
class Difference;


template< class EntityImp, class MinuendType, class SubtrahendType >
class DifferenceTraits
{
public:
  typedef typename MinuendType::template LocalFunction< EntityImp >::Type  LocalMinuendType;
  typedef typename SubtrahendType::template LocalFunction< EntityImp >::Type LocalSubtrahendType;
  typedef typename LocalMinuendType::DomainFieldType DomainFieldType;
  static_assert(std::is_same< DomainFieldType, typename LocalSubtrahendType::DomainFieldType >::value,
                "DomainFieldType of LocalMinuendType and LocalSubtrahendType do not match!");
  static const unsigned int dimDomain = LocalMinuendType::dimDomain;
  static_assert(dimDomain == LocalSubtrahendType::dimDomain,
                "dimDomain of LocalMinuendType and LocalSubtrahendType do not match!");
  typedef typename LocalMinuendType::RangeFieldType RangeFieldType;
  static_assert(std::is_same< RangeFieldType, typename LocalSubtrahendType::RangeFieldType >::value,
                "RangeFieldType of LocalMinuendType and LocalSubtrahendType do not match!");
  static const unsigned int dimRangeRows = LocalMinuendType::dimRangeRows;
  static_assert(dimRangeRows == LocalSubtrahendType::dimRangeRows,
                "dimRangeRows of LocalMinuendType and LocalSubtrahendType do not match!");
  static const unsigned int dimRangeCols = LocalMinuendType::dimRangeCols;
  static_assert(dimRangeCols == LocalSubtrahendType::dimRangeCols,
                "dimRangeCols of LocalMinuendType and LocalSubtrahendType do not match!");
  typedef EntityImp EntityType;
  typedef Difference< EntityType, MinuendType, SubtrahendType > derived_type;
  static_assert(std::is_base_of<  Stuff::LocalFunctionInterface<  typename LocalMinuendType::Traits,
                                                                  DomainFieldType, dimDomain,
                                                                  RangeFieldType, dimRangeRows, dimRangeCols >,
                                  LocalMinuendType >::value,
                "LocalMinuendType has to be derived from Stuff::LocalFunctionInterface!");
  static_assert(std::is_base_of<  Stuff::LocalFunctionInterface<  typename LocalSubtrahendType::Traits,
                                                                  DomainFieldType, dimDomain,
                                                                  RangeFieldType, dimRangeRows, dimRangeCols >,
                                  LocalSubtrahendType >::value,
                "LocalSubtrahendType has to be derived from Stuff::LocalFunctionInterface!");
};


template< class EntityImp, class MinuendType, class SubtrahendType  >
class Difference
  : public Stuff::LocalFunctionInterface< DifferenceTraits< EntityImp, MinuendType, SubtrahendType >,
                                          typename DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::DomainFieldType,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimDomain,
                                          typename DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::RangeFieldType,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimRangeRows,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimRangeCols >
{
  typedef Stuff::LocalFunctionInterface<  DifferenceTraits< EntityImp, MinuendType, SubtrahendType >,
                                          typename DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::DomainFieldType,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimDomain,
                                          typename DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::RangeFieldType,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimRangeRows,
                                          DifferenceTraits< EntityImp, MinuendType, SubtrahendType >::dimRangeCols > BaseType;
  typedef Difference< EntityImp, MinuendType, SubtrahendType > ThisType;
public:
  typedef DifferenceTraits< EntityImp, MinuendType, SubtrahendType >  Traits;
private:
  typedef typename Traits::LocalMinuendType LocalMinuendType;
  typedef typename Traits::LocalSubtrahendType LocalSubtrahendType;
public:
  typedef typename Traits::EntityType                                 EntityType;

  typedef typename BaseType::DomainType         DomainType;
  typedef typename BaseType::RangeType          RangeType;
  typedef typename BaseType::JacobianRangeType  JacobianRangeType;

  Difference(const EntityType& entity, const MinuendType& minuend, const SubtrahendType& subtrahend)
    : entity_(entity)
    , minuend_(minuend)
    , subtrahend_(subtrahend)
    , local_minuend_(new LocalMinuendType(minuend_.localFunction(entity_)))
    , local_subtrahend_(new LocalSubtrahendType(subtrahend_.localFunction(entity_)))
    , tmp_value_(0)
    , tmp_jacobian_value_(0)
  {}

  Difference(const ThisType& other)
    : entity_(other.entity_)
    , minuend_(other.minuend_)
    , subtrahend_(other.subtrahend_)
    , local_minuend_(new LocalMinuendType(minuend_.localFunction(entity_)))
    , local_subtrahend_(new LocalSubtrahendType(subtrahend_.localFunction(entity_)))
    , tmp_value_(0)
    , tmp_jacobian_value_(0)
  {}

  ~Difference()
  {
    delete local_minuend_;
    delete local_subtrahend_;
  }

private:
  ThisType& operator=(const ThisType& other);

public:
  const EntityType& entity() const
  {
    return entity_;
  }

  virtual int order() const
  {
    if ((local_minuend_->order() < 0) || (local_subtrahend_->order() < 0))
      return -1;
    else
      return std::max(local_minuend_->order(), local_subtrahend_->order());
  }

  void evaluate(const DomainType& xx, RangeType& ret) const
  {
    local_minuend_->evaluate(xx, ret);
    local_subtrahend_->evaluate(xx, tmp_value_);
    ret -= tmp_value_;
  }

  void jacobian(const DomainType& xx, JacobianRangeType& ret) const
  {
    local_minuend_->jacobian(xx, ret);
    local_subtrahend_->jacobian(xx, tmp_jacobian_value_);
    ret -= tmp_jacobian_value_;
  }

private:
  const EntityType& entity_;
  const MinuendType& minuend_;
  const SubtrahendType& subtrahend_;
  const LocalMinuendType* local_minuend_;
  const LocalSubtrahendType* local_subtrahend_;
  mutable RangeType tmp_value_;
  mutable JacobianRangeType tmp_jacobian_value_;
};


} // namespace LocalFunction

namespace DiscreteFunction {


template< class MinuendType, class SubtrahendType >
class Difference
  : public Stuff::LocalizableFunction
{
  static_assert(std::is_base_of< Stuff::LocalizableFunction, MinuendType >::value,
                "MinuendType has to be derived from Stuff::LocalizableFunction!");
  static_assert(std::is_base_of< Stuff::LocalizableFunction, SubtrahendType >::value,
                "SubtrahendType has to be derived from Stuff::LocalizableFunction!");
public:
  template< class EntityType >
  class LocalFunction
  {
  public:
    typedef GDT::LocalFunction::Difference< EntityType, MinuendType, SubtrahendType > Type;
  };

  Difference(const MinuendType& minuend, const SubtrahendType& subtrahend)
    : minuend_(minuend)
    , subtrahend_(subtrahend)
  {}

  template< class EntityType >
  typename LocalFunction< EntityType >::Type localFunction(const EntityType& entity) const
  {
    return typename LocalFunction< EntityType >::Type(entity, minuend_, subtrahend_);
  }

private:
  const MinuendType& minuend_;
  const SubtrahendType& subtrahend_;
}; // class Difference


} // namespace DiscreteFunction
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DIFFERENCE_HH
