// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_INTERFACES_HH
#define DUNE_GDT_PRODUCTS_INTERFACES_HH

#include <type_traits>

#include <dune/stuff/common/crtp.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Interface for all products (that are neither localizable or assemblable).
 * \todo  Add more documentation.
 */
template <class Traits>
class ProductInterface : public Stuff::CRTPInterface<ProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FieldType FieldType;

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp().grid_view());
    return this->as_imp().grid_view();
  }

  template <class RangeType, class SourceType>
  FieldType apply2(const RangeType& range, const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp().apply2(range, source));
    return this->as_imp().apply2(range, source);
  }

  template <class RangeType>
  FieldType induced_norm(const RangeType& range) const
  {
    return std::sqrt(apply2(range, range));
  }
}; // class ProductInterface


/**
 * \brief Interface for all localizable products.
 * \todo  Add more documentatin, especially about the notion of localibale.
 */
template <class Traits>
class LocalizableProductInterface : public Stuff::CRTPInterface<LocalizableProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::SourceType SourceType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

private:
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, SourceType>::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, RangeType>::value,
                "RangeType has to be derived from Stuff::IsLocalizableFunction!");
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
    CHECK_CRTP(this->as_imp().grid_view());
    return this->as_imp().grid_view();
  }

  const RangeType& range() const
  {
    CHECK_CRTP(this->as_imp().range());
    return this->as_imp().range();
  }

  const SourceType& source() const
  {
    CHECK_CRTP(this->as_imp().source());
    return this->as_imp().source();
  }

  FieldType apply2()
  {
    CHECK_CRTP(this->as_imp().apply2());
    return this->as_imp().apply2();
  }
}; // class LocalizableProductInterface


/**
 * \brief Interface for all assemblable products.
 * \todo  Add more documentatin, especially about the notion of assemblebla.
 */
template <class Traits>
class AssemblableProductInterface : public Stuff::CRTPInterface<AssemblableProductInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::MatrixType MatrixType;

  typedef typename MatrixType::ScalarType FieldType;
  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const size_t dimDomain = GridViewType::dimension;

  typedef Stuff::LA::SparsityPatternDefault PatternType;

private:
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceType::Traits, SourceSpaceType::dimDomain,
                                               SourceSpaceType::dimRange, SourceSpaceType::dimRangeCols>,
                                SourceSpaceType>::value,
                "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceType::Traits, RangeSpaceType::dimDomain,
                                               RangeSpaceType::dimRange, RangeSpaceType::dimRangeCols>,
                                RangeSpaceType>::value,
                "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same<typename RangeSpaceType::GridViewType, GridViewType>::value,
                "The GridViewType of RangeSpaceType and GridViewType have to match!");
  static_assert(std::is_same<typename SourceSpaceType::GridViewType, GridViewType>::value,
                "The GridViewType of SourceSpaceType and GridViewType have to match!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixType::Traits, FieldType>, MatrixType>::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");

public:
  static PatternType pattern(const RangeSpaceType& range_space)
  {
    return pattern(range_space, range_space);
  }

  static PatternType pattern(const RangeSpaceType& range_space, const SourceSpaceType& source_space)
  {
    return pattern(range_space, source_space, range_space.grid_view());
  }

  static PatternType pattern(const RangeSpaceType& range_space, const GridViewType& grid_view)
  {
    return pattern(range_space, range_space, grid_view);
  }

  static PatternType pattern(const RangeSpaceType& range_space, const SourceSpaceType& source_space,
                             const GridViewType& grid_view)
  {
    return derived_type::pattern(range_space, source_space, grid_view);
  }

  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp().grid_view());
    return this->as_imp().grid_view();
  }

  const RangeSpaceType& range_space() const
  {
    CHECK_CRTP(this->as_imp().range_space());
    return this->as_imp().range_space();
  }

  const SourceSpaceType& source_space() const
  {
    CHECK_CRTP(this->as_imp().source_space());
    return this->as_imp().source_space();
  }

  void assemble()
  {
    CHECK_AND_CALL_CRTP(this->as_imp().assemble());
  }

  MatrixType& matrix()
  {
    CHECK_CRTP(this->as_imp().matrix());
    return this->as_imp().matrix();
  }

  const MatrixType& matrix() const
  {
    CHECK_CRTP(this->as_imp().matrix());
    return this->as_imp().matrix();
  }

  template <class R, class S>
  FieldType apply2(const Stuff::LA::VectorInterface<R, FieldType>& range,
                   const Stuff::LA::VectorInterface<S, FieldType>& source)
  {
    assemble();
    assert(range.size() == matrix().rows());
    assert(source.size() == matrix().cols());
    auto tmp = range.copy();
    matrix().mv(source.as_imp(source), tmp);
    return range.dot(tmp);
  } // ... apply2(...)

  template <class R, class S>
  FieldType apply2(const ConstDiscreteFunction<RangeSpaceType, R>& range,
                   const ConstDiscreteFunction<SourceSpaceType, S>& source)
  {
    return apply2(range.vector(), source.vector());
  }

  template <class R>
  FieldType induced_norm(const Stuff::LA::VectorInterface<R, FieldType>& range)
  {
    return std::sqrt(apply2(range, range));
  }

  template <class R>
  FieldType induced_norm(const ConstDiscreteFunction<RangeSpaceType, R>& range)
  {
    return std::sqrt(apply2(range, range));
  }
}; // class AssemblableProductInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_INTERFACES_HH
