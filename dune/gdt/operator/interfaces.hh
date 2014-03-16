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

#include <dune/gdt/space/interface.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/assembler/gridwalker.hh>
#include <dune/gdt/la/containerfactory.hh>
#include <dune/gdt/assembler/local/codim0.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class ProductInterface
  : protected Stuff::CRTPInterface< ProductInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::FieldType    FieldType;

  const GridPartType& grid_part() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_part());
    return this->as_imp(*this).grid_part();
  }

  template< class RangeType, class SourceType >
  FieldType apply2(const RangeType& range, const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp(*this).apply2(range, source));
    return this->as_imp(*this).apply2(range, source);
  }
}; // class ProductInterface


template< class Traits >
class LocalizableProductInterface
  : protected Stuff::CRTPInterface< LocalizableProductInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::FieldType    FieldType;

  typedef typename GridPartType::template Codim< 0 >::EntityType  EntityType;
  typedef typename GridPartType::ctype                            DomainFieldType;
  static const unsigned int                                       dimDomain = GridPartType::dimension;

private:
  static_assert(std::is_base_of< Stuff::IsLocalizableFunction, SourceType >::value,
                "SourceType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_base_of< Stuff::IsLocalizableFunction, RangeType >::value,
                "RangeType has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same< typename SourceType::EntityType, EntityType >::value,
                "The EntityType of SourceType and GridPartType have to match!");
  static_assert(std::is_same< typename RangeType::EntityType, EntityType >::value,
                "The EntityType of RangeType and GridPartType have to match!");
  static_assert(std::is_same< typename SourceType::DomainFieldType, DomainFieldType >::value,
                "The DomainFieldType of SourceType and GridPartType have to match!");
  static_assert(std::is_same< typename RangeType::DomainFieldType, DomainFieldType >::value,
                "The DomainFieldType of RangeType and GridPartType have to match!");
  static_assert(SourceType::dimDomain == dimDomain, "The dimDomain of SourceType and GridPartType have to match!");
  static_assert(RangeType::dimDomain == dimDomain, "The dimDomain of RangeType and GridPartType have to match!");

public:
  const GridPartType& grid_part() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_part());
    return this->as_imp(*this).grid_part();
  }

  const RangeType& range() const
  {
    CHECK_CRTP(this->as_imp(*this).range());
    return this->as_imp(*this).range();
  }

  const SourceType& source() const
  {
    CHECK_CRTP(this->as_imp(*this).source());
    return this->as_imp(*this).source();
  }

  FieldType apply2()
  {
    CHECK_CRTP(this->as_imp(*this).apply2());
    return this->as_imp(*this).apply2();
  }
}; // class LocalizableProductInterface


// forward, to be used in the traits
template< class ImpTraits >
class LocalizableProductBase;

template< class ImpTraits >
class LocalizableProductBaseTraits
{
public:
  typedef typename ImpTraits::derived_type  derived_type;
  typedef typename ImpTraits::GridPartType  GridPartType;
  typedef typename ImpTraits::RangeType     RangeType;
  typedef typename ImpTraits::SourceType    SourceType;
  typedef typename ImpTraits::FieldType     FieldType;
}; // class LocalizableProductBaseTraits


template< class ImpTraits >
class LocalizableProductBase
  : public LocalizableProductInterface< LocalizableProductBaseTraits< ImpTraits > >
  , public LocalizableProductInterface< ImpTraits >
  , public Functor::Codim0< typename ImpTraits::GridPartType >
{
  typedef LocalizableProductInterface< LocalizableProductBaseTraits< ImpTraits > > InterfaceType;
  typedef LocalizableProductBaseTraits< ImpTraits > Traits;
public:
  typedef typename Traits::GridPartType GridPartType;
  typedef typename Traits::RangeType    RangeType;
  typedef typename Traits::SourceType   SourceType;
  typedef typename Traits::FieldType    FieldType;
private:
  typedef typename ImpTraits::LocalOperatorType LocalOperatorType;

public:
  using typename InterfaceType::EntityType;

public:
  LocalizableProductBase(const GridPartType& grid_part, const RangeType& range, const SourceType& source)
    : grid_part_(grid_part)
    , range_(range)
    , source_(source)
    , prepared_(false)
    , tmp_local_operator_result_(1, 1, 0)
    , tmp_matrices_()
    , result_(0)
  {}

  virtual ~LocalizableProductBase() {}

  const GridPartType& grid_part() const
  {
    return grid_part_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  const SourceType& source() const
  {
    return source_;
  }

private:
  virtual const LocalOperatorType& local_operator() const = 0;

public:
  virtual void prepare() DS_OVERRIDE
  {
    if (!prepared_) {
      tmp_matrices_ = std::vector< DynamicMatrix< FieldType > >(local_operator().numTmpObjectsRequired(),
                                                                DynamicMatrix< FieldType >(1, 1, 0));
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE
  {
    assert(prepared_);
    // get the local functions
    const auto local_source_ptr = this->source().local_function(entity);
    const auto local_range_ptr = this->range().local_function(entity);
    // apply local operator
    local_operator().apply(*local_source_ptr, *local_range_ptr, tmp_local_operator_result_, tmp_matrices_);
    assert(tmp_local_operator_result_.rows() == 1);
    assert(tmp_local_operator_result_.cols() == 1);
    result_ += tmp_local_operator_result_[0][0];
  } // ... apply_local(...)

  FieldType apply2()
  {
    result_ *= FieldType(0);
    GridWalker< GridPartType > grid_walker(grid_part_);
    grid_walker.add(*this);
    grid_walker.walk();
    return result_;
  }

private:
  const GridPartType& grid_part_;
  const RangeType& range_;
  const SourceType& source_;
  bool prepared_;
  DynamicMatrix< FieldType > tmp_local_operator_result_;
  std::vector< DynamicMatrix< FieldType > > tmp_matrices_;
  FieldType result_;
}; // class LocalizableProductBase


template< class Traits >
class AssemblableProductInterface
  : protected Stuff::CRTPInterface< AssemblableProductInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type     derived_type;
  typedef typename Traits::GridPartType     GridPartType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::MatrixType       MatrixType;

  typedef typename MatrixType::ScalarType                         FieldType;
  typedef typename GridPartType::template Codim< 0 >::EntityType  EntityType;
  typedef typename GridPartType::ctype                            DomainFieldType;
  static const unsigned int                                       dimDomain = GridPartType::dimension;

private:
  static_assert(std::is_base_of< SpaceInterface< typename RangeSpaceType::Traits >, RangeSpaceType >::value,
                "RangeSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of< SpaceInterface< typename SourceSpaceType::Traits >, SourceSpaceType >::value,
                "SourceSpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_same< typename RangeSpaceType::GridPartType, GridPartType >::value,
                "The GridPartType of RangeSpaceType and GridPartType have to match!");
  static_assert(std::is_same< typename SourceSpaceType::GridPartType, GridPartType >::value,
                "The GridPartType of SourceSpaceType and GridPartType have to match!");
  static_assert(std::is_base_of< Stuff::LA::MatrixInterface< typename MatrixType::Traits >, MatrixType >::value,
                "MatrixType has to be derived from Stuff::LA::MatrixInterface!");

public:
  const GridPartType& grid_part() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_part());
    return this->as_imp(*this).grid_part();
  }

  const RangeSpaceType& range_space() const
  {
    CHECK_CRTP(this->as_imp(*this).range_space());
    return this->as_imp(*this).range_space();
  }

  const SourceSpaceType& source_space() const
  {
    CHECK_CRTP(this->as_imp(*this).source_space());
    return this->as_imp(*this).source_space();
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

  template< class R, class S >
  FieldType apply2(const Stuff::LA::VectorInterface< R >& range,
                   const Stuff::LA::VectorInterface< S >& source)
  {
    CHECK_CRTP(this->as_imp(*this).apply2(range, source));
    return this->as_imp(*this).apply2(range, source);
  }

  template< class R, class S >
  FieldType apply2(const ConstDiscreteFunction< RangeSpaceType, R >& range,
                   const ConstDiscreteFunction< SourceSpaceType, S >& source)
  {
    apply2(range.vector(), source.vector());
  }
}; // class AssemblableProductInterface


// forward, to be used in the traits
template< class ImpTraits >
class AssemblableProductBase;


template< class ImpTraits >
class AssemblableProductBaseTraits
{
public:
  typedef typename ImpTraits::derived_type    derived_type;
  typedef typename ImpTraits::GridPartType    GridPartType;
  typedef typename ImpTraits::RangeSpaceType  RangeSpaceType;
  typedef typename ImpTraits::SourceSpaceType SourceSpaceType;
  typedef typename ImpTraits::MatrixType      MatrixType;
}; // class AssemblableProductBaseTraits


template< class ImpTraits >
class AssemblableProductBase
    : public AssemblableProductInterface< AssemblableProductBaseTraits< ImpTraits > >
    , public AssemblableProductInterface< ImpTraits >
    , public Functor::Codim0< typename ImpTraits::GridPartType >
{
  typedef AssemblableProductInterface< AssemblableProductBaseTraits< ImpTraits > > InterfaceType;
  typedef AssemblableProductBaseTraits< ImpTraits > Traits;
public:
  typedef typename Traits::GridPartType     GridPartType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::MatrixType       MatrixType;

  typedef typename MatrixType::ScalarType   FieldType;
private:
  typedef typename ImpTraits::LocalOperatorType LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix< LocalOperatorType > LocalAssemblerType;

public:
  using typename InterfaceType::EntityType;

  AssemblableProductBase(const GridPartType& grid_part,
                         const RangeSpaceType& range_space,
                         const SourceSpaceType& source_space)
    : grid_part_(grid_part)
    , range_space_(range_space)
    , source_space_(source_space)
    , matrix_(LA::ContainerFactory< MatrixType >::create(range_space, source_space))
    , local_assembler_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {}

  const GridPartType& grid_part() const
  {
    return grid_part_;
  }

  const RangeSpaceType& range_space() const
  {
    return range_space_;
  }

  const SourceSpaceType& source_space() const
  {
    return source_space_;
  }

  MatrixType& matrix()
  {
    return matrix_;
  }

  const MatrixType& matrix() const
  {
    return matrix_;
  }

private:
  virtual const LocalOperatorType& local_operator() const = 0;

public:
  virtual void prepare()
  {
    if (!assembled_ && !prepared_) {
      local_assembler_ = std::unique_ptr< LocalAssemblerType >(new LocalAssemblerType(local_operator()));
      const auto num_tmp_objects_required = local_assembler_->numTmpObjectsRequired();
      const size_t max_local_size = std::max(range_space_.mapper().maxNumDofs(), source_space_.mapper().maxNumDofs());
      tmp_local_matrices_
          = { std::vector< DynamicMatrix< FieldType > >(num_tmp_objects_required[0],
                                                        DynamicMatrix< FieldType >(max_local_size,
                                                                                   max_local_size))
            , std::vector< DynamicMatrix< FieldType > >(num_tmp_objects_required[1],
                                                        DynamicMatrix< FieldType >(max_local_size,
                                                                                   max_local_size))
            };
      tmp_indices_container_ = { DynamicVector< size_t >(max_local_size, 0)
                               , DynamicVector< size_t >(max_local_size, 0)
                               };
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity)
  {
    assert(prepared_);
    local_assembler_->assembleLocal(range_space_, source_space_,
                                    entity,
                                    matrix_,
                                    tmp_local_matrices_, tmp_indices_container_);
  } // ... apply_local(...)

  void assemble()
  {
    if (!assembled_) {
      GridWalker< GridPartType > grid_walker(grid_part_);
      grid_walker.add(*this);
      grid_walker.walk();
      assembled_ = true;
    }
  } // ... assemble()

  template< class R, class S >
  FieldType apply2(const Stuff::LA::VectorInterface< R >& range, const Stuff::LA::VectorInterface< S >& source)
  {
    assemble();
    auto tmp = range.copy();
    matrix_.mv(source, tmp);
    return range.dot(tmp);
  } // ... apply2(...)

private:
  const GridPartType& grid_part_;
  const RangeSpaceType& range_space_;
  const SourceSpaceType& source_space_;
  MatrixType matrix_;
  std::shared_ptr< LocalAssemblerType > local_assembler_;
  bool prepared_;
  bool assembled_;
  std::vector< std::vector< DynamicMatrix< FieldType > > > tmp_local_matrices_;
  std::vector< DynamicVector< size_t > > tmp_indices_container_;
}; // class AssemblableProductBase


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATOR_INTERFACES_HH
