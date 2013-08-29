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

template< class MinuendType, class SubrahendImp, class D, int d, class R, int rR, int rC = 1 >
class Difference;

template< class MinuendType, class SubrahendImp, class D, int d, class R, int rR, int rC = 1 >
class DifferenceTraits;


template< class MinuendType, class SubtrahendType, class D, int d, class R, int r >
class DifferenceTraits< MinuendType, SubtrahendType, D, d, R, r, 1 >
{
  static_assert(std::is_base_of< Stuff::LocalFunctionInterface< typename MinuendType::Traits, D, d, R, r, 1 >,
                                                                MinuendType >::value,
                "MinuendType has to be derived from Stuff::LocalFunctionInterface< ..., D, d, R, r, 1 >");
  static_assert(std::is_base_of< Stuff::LocalFunctionInterface< typename SubtrahendType::Traits, D, d, R, r, 1 >,
                                                                SubtrahendType >::value,
                "SubtrahendType has to be derived from Stuff::LocalFunctionInterface< ..., D, d, R, r, 1 >");
  static_assert(std::is_same< typename MinuendType::EntityType, typename SubtrahendType::EntityType >::value,
                "EntityType of MinuendType and SubtrahendType do not match!");
public:
  typedef Difference< MinuendType, SubtrahendType, D, d, R, r, 1 > derived_type;
  typedef typename MinuendType::EntityType EntityType;
};


template< class MinuendType, class SubtrahendType, class D, int d, class R, int r >
class Difference< MinuendType, SubtrahendType, D, d, R, r, 1 >
  : public Stuff::LocalFunctionInterface< DifferenceTraits< MinuendType, SubtrahendType, D, d, R, r, 1 >, D, d, R, r, 1 >
{
public:
  typedef DifferenceTraits< MinuendType, SubtrahendType, D, d, R, r, 1 > Traits;
  typedef typename Traits::EntityType EntityType;

  typedef D                                               DomainFieldType;
  static const unsigned int                               dimDomain = d;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;

  typedef R                                             RangeFieldType;
  static const unsigned int                             dimRange = r;
  static const unsigned int                             dimRangeRows = dimRange;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;

  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  Difference(const MinuendType& minuend, const SubtrahendType& subtrahend)
    : minuend_(minuend)
    , subtrahend_(subtrahend)
    , tmp_value_(0)
    , tmp_jacobian_value_(0)
  {
    assert(minuend_.entity() == subtrahend_.entity());
  }

  const EntityType& entity() const
  {
    return minuend_.entity();
  }

  virtual int order() const
  {
    if ((minuend_.order() < 0) || (subtrahend_.order() < 0))
      return -1;
    else
      return std::max(minuend_.order(), subtrahend_.order());
  }

  void evaluate(const DomainType& xx, RangeType& ret) const
  {
    minuend_.evaluate(xx, ret);
    subtrahend_.evaluate(xx, tmp_value_);
    ret -= tmp_value_;
  }

  void jacobian(const DomainType& xx, JacobianRangeType& ret) const
  {
    minuend_.evaluate(xx, ret);
    subtrahend_.evaluate(xx, tmp_jacobian_value_);
    ret -= tmp_jacobian_value_;
  }

private:
  const MinuendType& minuend_;
  const SubtrahendType& subtrahend_;
  mutable RangeType tmp_value_;
  mutable JacobianRangeType tmp_jacobian_value_;
};


} // namespace LocalFunction

namespace DiscreteFunction {


template< class MinuendType, class SubtrahendType >
class Difference
  : public Stuff::LocalizableFunction
{
public:
  template< class EntityType >
  class LocalFunction
  {
    typedef typename MinuendType::template LocalFunction< EntityType >     LocalMinuendType;
    typedef typename SubtrahendType::template LocalFunction< EntityType >  LocalSubtrahendType;
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
  public:
    typedef GDT::LocalFunction::Difference< LocalMinuendType,
                                            LocalSubtrahendType,
                                            DomainFieldType, dimDomain,
                                            RangeFieldType, dimRangeRows, dimRangeCols > Type;
  };

  Difference(const MinuendType& minuend, const SubtrahendType& subtrahend)
    : minuend_(minuend)
    , subtrahend_(subtrahend)
  {}

  template< class EntityType >
  typename LocalFunction< EntityType >::Type localFunction(const EntityType& entity) const
  {
    return typename LocalFunction< EntityType >::Type(minuend_.localFunction(entity),
                                                      subtrahend_.localFunction(entity));
  }

private:
  const MinuendType& minuend_;
  const SubtrahendType& subtrahend_;
}; // class Difference


} // namespace DiscreteFunction
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_DIFFERENCE_HH
