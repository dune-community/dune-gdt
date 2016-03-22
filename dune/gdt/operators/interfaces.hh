// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_INTERFACES_HH
#define DUNE_GDT_OPERATORS_INTERFACES_HH

#include <type_traits>

#include <dune/common/deprecated.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

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


/**
 * \note The methods apply() and apply2() do not support threaded assembly atm, else we would have to add a use_tbb
 *       switch (to a lot of methods). Either this gets merged with the new Parameter or we need a different paradigm.
 */
template< class Traits >
class OperatorInterface
  : public Stuff::CRTPInterface< OperatorInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::FieldType    FieldType;
//  typedef typename Traits::JacobianType JacobianType;

  /// \name Methods that have to be implemented by any derived class
  /// \{

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range) const
  {
    CHECK_CRTP(this->as_imp().apply(source, range));
  }

  template< class RangeType, class SourceType >
  FieldType apply2(const RangeType& range, const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp().apply2(range, source));
    return this->as_imp().apply2(range, source);
  }

//  template< class SourceType >
//  JacobianType jacobian(const SourceType& source) const
//  {
//    CHECK_CRTP(this->as_imp().jacobian(source));
//    return this->as_imp().jacobian(source);
//  }

  template< class RangeType, class SourceType >
  void apply_inverse(const RangeType& range, SourceType& source, const Stuff::Common::Configuration& opts) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().apply_inverse(range, source, opts));
  }

  std::vector< std::string > invert_options() const
  {
    CHECK_CRTP(this->as_imp().invert_options());
    return this->as_imp().invert_options();
  }

  Stuff::Common::Configuration invert_options(const std::string& type) const
  {
    CHECK_CRTP(this->as_imp().invert_options(type));
    return this->as_imp().invert_options(type);
  }

  /// \}
  /// \name Provided by the interface for convenience
  /// \{

  template< class RangeType, class SourceType >
  void apply_inverse(const RangeType& range, SourceType& source, const std::string& type) const
  {
    apply_inverse(range, source, invert_options(type));
  }

  template< class RangeType, class SourceType >
  void apply_inverse(const RangeType& range, SourceType& source) const
  {
    auto type = invert_options();
    apply_inverse(range, source, type.size() > 0 ? type[0] : "");
  }

  template< class RangeType >
  FieldType induced_norm(const RangeType& range) const
  {
    return std::sqrt(apply2(range, range));
  }

  /// \}
}; // class OperatorInterface


template< class Traits >
class
    DUNE_DEPRECATED_MSG("Derive from LocalizableOperatorDefault instead (04.07.2015)!")
      LocalizableOperatorInterface
  : public Stuff::CRTPInterface< LocalizableOperatorInterface< Traits >, Traits >
{
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridViewType::template Codim< 0 >::Entity EntityType;
  typedef typename GridViewType::ctype                       DomainFieldType;
  static const size_t                                        dimDomain = GridViewType::dimension;

private:
  static_assert(Stuff::is_localizable_function< SourceType >::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(is_discrete_function< RangeType >::value, "RangeType has to be derived from DiscreteFunction!");
  static_assert(std::is_same< typename SourceType::EntityType, EntityType >::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::EntityType, EntityType >::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same< typename SourceType::DomainFieldType, DomainFieldType >::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::DomainFieldType, DomainFieldType >::value,
                "The DomainFieldType of RangeType and GridViewType have to match!");
  static_assert(SourceType::dimDomain == dimDomain, "The dimDomain of SourceType and GridViewType have to match!");
  static_assert(RangeType::dimDomain == dimDomain,  "The dimDomain of RangeType and GridViewType have to match!");

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


template< class Traits >
class
    DUNE_DEPRECATED_MSG("Derive from MatrixOperatorDefault instead (04.07.2015)!")
      AssemblableOperatorInterface
  : public AssemblableProductInterface< Traits >
{
  typedef AssemblableProductInterface< Traits > BaseType;
public:
  typedef typename BaseType::derived_type derived_type;

  using typename BaseType::FieldType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::DomainFieldType;

  template< class S, class R >
  void apply(const Stuff::LA::VectorInterface< S, DomainFieldType >& source,
             Stuff::LA::VectorInterface< R, FieldType >& range)
  {
    CHECK_CRTP(this->as_imp().apply(source.as_imp(), range.as_imp()));
    return this->as_imp().apply(source.as_imp(), range.as_imp());
  }

  template< class S, class R >
  void apply(const ConstDiscreteFunction< SourceSpaceType, S >& source,
             ConstDiscreteFunction< RangeSpaceType, R >& range)
  {
    apply(source.vector(), range.vector());
  }

  static std::vector< std::string > invert_options()
  {
    return derived_type::invert_options();
  }

  static Stuff::Common::Configuration invert_options(const std::string& type)
  {
    return derived_type::invert_options(type);
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R, FieldType >& range,
                     const Stuff::LA::VectorInterface< S, FieldType >& source)
  {
    apply_inverse(range, source, invert_options()[0]);
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R, FieldType >& range,
                     Stuff::LA::VectorInterface< S, FieldType >& source,
                     const std::string& opt)
  {
    apply_inverse(range, source, invert_options(opt));
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R, FieldType >& range,
                     Stuff::LA::VectorInterface< S, FieldType >& source,
                     const Stuff::Common::Configuration& opts)
  {
    CHECK_CRTP(this->as_imp(*this).apply_inverse(range.as_imp(range), source.as_imp(source), opts));
    return this->as_imp(*this).apply_inverse(range.as_imp(range), source.as_imp(source), opts);
  }

  template< class R, class S >
  void apply_inverse(const ConstDiscreteFunction< SourceSpaceType, R >& range,
                     ConstDiscreteFunction< RangeSpaceType, S >& source)
  {
    apply_inverse(range.vector(), source.vector());
  }

  template< class R, class S >
  void apply_inverse(const ConstDiscreteFunction< SourceSpaceType, R >& range,
                     ConstDiscreteFunction< RangeSpaceType, S >& source,
                     const std::string& opt)
  {
    apply_inverse(range.vector(), source.vector(), opt);
  }

  template< class R, class S >
  void apply_inverse(const ConstDiscreteFunction< SourceSpaceType, R >& range,
                     ConstDiscreteFunction< RangeSpaceType, S >& source,
                     const Stuff::Common::Configuration& opts)
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }

  template< class R, class S, class P >
  FieldType apply2(const Stuff::LA::VectorInterface< R, FieldType >& range,
                   const Stuff::LA::VectorInterface< S, FieldType >& source,
                   const AssemblableProductInterface< P >& product)
  {
    auto tmp = range.copy();
    apply(source, tmp);
    return product.apply2(tmp, source);
  }

  template< class R, class S >
  FieldType apply2(const Stuff::LA::VectorInterface< R, FieldType >& range,
                   const Stuff::LA::VectorInterface< S, FieldType >& source)
  {
    this->assemble();
    auto tmp = range.copy();
    this->matrix().mv(source.as_imp(), tmp);
    return range.dot(tmp);
  }

  template< class S, class R >
  FieldType apply2(const ConstDiscreteFunction< SourceSpaceType, S >& source,
                   const ConstDiscreteFunction< RangeSpaceType, R >& range)
  {
    return apply2(range.vector(), source.vector());
  }
}; // class AssemblableOperatorInterface


namespace internal {


template< class Tt >
struct is_operator_helper
{
  DSC_has_typedef_initialize_once(Traits)

  static const bool is_candidate = DSC_has_typedef(Traits)< Tt >::value;
};


} // namespace internal


template< class T, bool candidate = internal::is_operator_helper< T >::is_candidate >
struct is_operator
  : public std::is_base_of< OperatorInterface< typename T::Traits >, T >
{};

template< class T >
struct is_operator< T, false >
  : public std::false_type
{};


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
