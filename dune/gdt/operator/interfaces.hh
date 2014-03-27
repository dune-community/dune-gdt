// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATOR_INTERFACES_HH
#define DUNE_GDT_OPERATOR_INTERFACES_HH

#include <type_traits>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/la/container/pattern.hh>
#include <dune/stuff/la/solver.hh>
#include <dune/stuff/common/configtree.hh>

#include <dune/gdt/space/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class OperatorInterface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::FieldType    FieldType;
  typedef typename Traits::GridViewType GridViewType;

  template< class SourceType, class RangeType >
  void apply(const SourceType& source, RangeType& range) const
  {
    CHECK_CRTP(this->as_imp(*this).apply(source, range));
    return this->as_imp(*this).apply(source, range);
  }
}; // class OperatorInterface


template< class Traits >
class LocalizableOperatorInterface
  : protected Stuff::CRTPInterface< LocalizableOperatorInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::RangeType    RangeType;

  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const unsigned int                                   dimDomain = GridViewType::dimension;

private:
  static_assert(std::is_base_of< Stuff::IsLocalizableFunction, SourceType >::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of< DiscreteFunction< typename RangeType::SpaceType, typename RangeType::VectorType >
                               , RangeType >::value,
                "RangeType has to be derived from DiscreteFunction!");
  static_assert(std::is_same< typename SourceType::EntityType, EntityType >::value,
                "The EntityType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::EntityType, EntityType >::value,
                "The EntityType of RangeType and GridViewType have to match!");
  static_assert(std::is_same< typename SourceType::DomainFieldType, DomainFieldType >::value,
                "The DomainFieldType of SourceType and GridViewType have to match!");
  static_assert(std::is_same< typename RangeType::DomainFieldType, DomainFieldType >::value,
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

  RangeType& range()
  {
    CHECK_CRTP(this->as_imp(*this).range());
    return this->as_imp(*this).range();
  }

  const RangeType& range() const
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
class AssemblableOperatorInterface
  : protected Stuff::CRTPInterface< AssemblableOperatorInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type     derived_type;
  typedef typename Traits::GridViewType     GridViewType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::MatrixType       MatrixType;

  typedef typename MatrixType::ScalarType                     FieldType;
  typedef typename GridViewType::template Codim< 0 >::Entity  EntityType;
  typedef typename GridViewType::ctype                        DomainFieldType;
  static const unsigned int                                   dimDomain = GridViewType::dimension;

private:
  static_assert(std::is_base_of< SpaceInterface< typename RangeSpaceType::Traits >, RangeSpaceType >::value,
                "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename SourceSpaceType::Traits >, SourceSpaceType >::value,
                "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same< typename RangeSpaceType::GridViewType, GridViewType >::value,
                "The GridViewType of RangeSpaceType and GridViewType have to match!");
  static_assert(std::is_same< typename SourceSpaceType::GridViewType, GridViewType >::value,
                "The GridViewType of SourceSpaceType and GridViewType have to match!");
  static_assert(std::is_base_of< Stuff::LA::MatrixInterface< typename MatrixType::Traits >, MatrixType >::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");

public:
  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space)
  {
    return pattern(range_space, range_space);
  }

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space)
  {
    return pattern(range_space, source_space, *(range_space.grid_view()));
  }

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return derived_type::pattern(range_space, source_space, grid_view);
  }

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  const SourceSpaceType& source_space() const
  {
    CHECK_CRTP(this->as_imp(*this).source_space());
    return this->as_imp(*this).source_space();
  }

  const RangeSpaceType& range_space() const
  {
    CHECK_CRTP(this->as_imp(*this).range_space());
    return this->as_imp(*this).range_space();
  }

  void assemble()
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).assemble());
  }

  MatrixType& matrix()
  {
    CHECK_CRTP(this->as_imp(*this).matrix());
    return this->as_imp(*this).matrix();
  }

  const MatrixType& matrix() const
  {
    CHECK_CRTP(this->as_imp(*this).matrix());
    return this->as_imp(*this).matrix();
  }

  template< class S, class R >
  void apply(const Stuff::LA::VectorInterface< S >& source, Stuff::LA::VectorInterface< R >& range)
  {
    CHECK_CRTP(this->as_imp(*this).apply(source, range));
    return this->as_imp(*this).apply(source, range);
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

  static Stuff::Common::ConfigTree invert_options(const std::string& type)
  {
    return derived_type::invert_options(type);
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R >& range, const Stuff::LA::VectorInterface< S >& source)
  {
    apply_inverse(range, source, invert_options()[0]);
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R >& range,
                     Stuff::LA::VectorInterface< S >& source,
                     const std::string& opt)
  {
    apply_inverse(range, source, invert_options(opt));
  }

  template< class R, class S >
  void apply_inverse(const Stuff::LA::VectorInterface< R >& range,
                     Stuff::LA::VectorInterface< S >& source,
                     const Stuff::Common::ConfigTree& opts)
  {
    CHECK_CRTP(this->as_imp(*this).apply_inverse(range, source, opts));
    return this->as_imp(*this).apply_inverse(range, source, opts);
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
                     const Stuff::Common::ConfigTree& opts)
  {
    apply_inverse(range.vector(), source.vector(), opts);
  }
}; // class AssemblableOperatorInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_INTERFACES_HH
