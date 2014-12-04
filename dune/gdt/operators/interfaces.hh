// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_INTERFACES_HH
#define DUNE_GDT_OPERATORS_INTERFACES_HH

#include <type_traits>

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/common/configuration.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

#include "../products/interfaces.hh"

namespace Dune {
namespace GDT {


template <class Traits>
class OperatorInterface : public ProductInterface<Traits>
{
  typedef ProductInterface<Traits> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::FieldType;

  template <class SourceType, class RangeType>
  void apply(const SourceType& source, RangeType& range) const
  {
    CHECK_CRTP(this->as_imp(*this).apply(source, range));
    return this->as_imp(*this).apply(source, range);
  }

  template <class S, class R>
  FieldType apply2(const Stuff::LA::VectorInterface<S, FieldType>& source,
                   const Stuff::LA::VectorInterface<R, FieldType>& range) const
  {
    auto tmp = range.copy();
    apply(source.as_imp(source), tmp);
    return range.dot(tmp);
  } // ... apply2(...)

  template <class SS, class SV, class RS, class RV>
  FieldType apply2(const ConstDiscreteFunction<SS, SV>& source, const ConstDiscreteFunction<RS, RV>& range) const
  {
    auto tmp_vector = range.copy();
    DiscreteFunction<RS, RV> tmp_function(range.space(), tmp_vector);
    apply(source, tmp_function);
    return range.vector().dot(tmp_vector);
  }
}; // class OperatorInterface


template <class Traits>
class LocalizableOperatorInterface : public Stuff::CRTPInterface<LocalizableOperatorInterface<Traits>, Traits>
{
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType SourceType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

private:
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, SourceType>::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<DiscreteFunction<typename RangeType::SpaceType, typename RangeType::VectorType>,
                                RangeType>::value,
                "RangeType has to be derived from DiscreteFunction!");
  static_assert(std::is_same<typename SourceType::EntityType, EntityType>::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::EntityType, EntityType>::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceType::DomainFieldType, DomainFieldType>::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same<typename RangeType::DomainFieldType, DomainFieldType>::value,
                "The DomainFieldType of RangeType and GridViewType have to match!");
  static_assert(SourceType::dimDomain == dimDomain, "The dimDomain of SourceType and GridViewType have to match!");
  static_assert(RangeType::dimDomain == dimDomain, "The dimDomain of RangeType and GridViewType have to match!");

public:
  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  const SourceType& source() const
  {
    CHECK_CRTP(this->as_imp(*this).source());
    return this->as_imp(*this).source();
  }

  const RangeType& range() const
  {
    CHECK_CRTP(this->as_imp(*this).range());
    return this->as_imp(*this).range();
  }

  RangeType& range()
  {
    CHECK_CRTP(this->as_imp(*this).range());
    return this->as_imp(*this).range();
  }

  void apply()
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).apply());
  }
}; // class LocalizableOperatorInterface


template <class Traits>
class AssemblableOperatorInterface : public AssemblableProductInterface<Traits>
{
  typedef AssemblableProductInterface<Traits> BaseType;

public:
  typedef typename BaseType::derived_type derived_type;

  using typename BaseType::GridViewType;
  using typename BaseType::FieldType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::MatrixType;

  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using BaseType::dimDomain;
  using typename BaseType::PatternType;

  template <class S, class R>
  void apply(const Stuff::LA::VectorInterface<S, DomainFieldType>& source,
             Stuff::LA::VectorInterface<R, FieldType>& range)
  {
    CHECK_CRTP(this->as_imp(*this).apply(source.as_imp(source), range.as_imp(range)));
    return this->as_imp(*this).apply(source.as_imp(source), range.as_imp(range));
  }

  template <class S, class R>
  void apply(const ConstDiscreteFunction<SourceSpaceType, S>& source, ConstDiscreteFunction<RangeSpaceType, R>& range)
  {
    apply(source.vector(), range.vector());
  }

  static std::vector<std::string> invert_options()
  {
    return derived_type::invert_options();
  }

  static Stuff::Common::Configuration invert_options(const std::string& type)
  {
    return derived_type::invert_options(type);
  }

  template <class R, class S>
  void apply_inverse(const Stuff::LA::VectorInterface<R, FieldType>& range,
                     const Stuff::LA::VectorInterface<S, FieldType>& source)
  {
    apply_inverse(range, source, invert_options()[0]);
  }

  template <class R, class S>
  void apply_inverse(const Stuff::LA::VectorInterface<R, FieldType>& range,
                     Stuff::LA::VectorInterface<S, FieldType>& source, const std::string& opt)
  {
    apply_inverse(range, source, invert_options(opt));
  }

  template <class R, class S>
  void apply_inverse(const Stuff::LA::VectorInterface<R, FieldType>& range,
                     Stuff::LA::VectorInterface<S, FieldType>& source, const Stuff::Common::Configuration& opts)
  {
    CHECK_CRTP(this->as_imp(*this).apply_inverse(range.as_imp(range), source.as_imp(source), opts));
    return this->as_imp(*this).apply_inverse(range.as_imp(range), source.as_imp(source), opts);
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<SourceSpaceType, R>& range,
                     ConstDiscreteFunction<RangeSpaceType, S>& source)
  {
    apply_inverse(range.vector(), source.vector());
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<SourceSpaceType, R>& range,
                     ConstDiscreteFunction<RangeSpaceType, S>& source, const std::string& opt)
  {
    apply_inverse(range.vector(), source.vector(), opt);
  }

  template <class R, class S>
  void apply_inverse(const ConstDiscreteFunction<SourceSpaceType, R>& range,
                     ConstDiscreteFunction<RangeSpaceType, S>& source, const Stuff::Common::Configuration& opts)
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }

  template <class R, class S, class P>
  FieldType apply2(const Stuff::LA::VectorInterface<R, FieldType>& range,
                   const Stuff::LA::VectorInterface<S, FieldType>& source,
                   const AssemblableProductInterface<P>& product)
  {
    auto tmp = range.copy();
    apply(source, tmp);
    return product.apply2(tmp, source);
  }

  template <class R, class S>
  FieldType apply2(const Stuff::LA::VectorInterface<R, FieldType>& range,
                   const Stuff::LA::VectorInterface<S, FieldType>& source)
  {
    this->assemble();
    auto tmp = range.copy();
    this->matrix().mv(source.as_imp(), tmp);
    return range.dot(tmp);
  }

  template <class S, class R>
  FieldType apply2(const ConstDiscreteFunction<SourceSpaceType, S>& source,
                   const ConstDiscreteFunction<RangeSpaceType, R>& range)
  {
    return apply2(range.vector(), source.vector());
  }
}; // class AssemblableOperatorInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
