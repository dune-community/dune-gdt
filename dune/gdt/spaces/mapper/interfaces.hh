// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_MAPPER_INTERFACES_HH
#define DUNE_GDT_SPACES_MAPPER_INTERFACES_HH

#include <dune/common/dynvector.hh>
#include <dune/grid/common/entity.hh>

#include <dune/xt/common/crtp.hh>
#include <dune/xt/common/type_traits.hh>

namespace Dune {
namespace GDT {


/**
 * \todo Add note that each implementation has to document which of globalIndices or mapToGlobal is more efficient!
 */
template <class Traits>
class MapperInterface : public XT::Common::CRTPInterface<MapperInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp(*this).backend());
    return this->as_imp(*this).backend();
  }

  size_t size() const
  {
    CHECK_CRTP(this->as_imp(*this).size());
    return this->as_imp(*this).size();
  }

  size_t maxNumDofs() const
  {
    CHECK_CRTP(this->as_imp(*this).maxNumDofs());
    return this->as_imp(*this).maxNumDofs();
  }

  template <int cd, class GridImp, template <int, int, class> class EntityImp>
  size_t numDofs(const Entity<cd, EntityType::dimension, GridImp, EntityImp>& entity) const
  {
    //    static_assert(std::is_same<Entity<cd, EntityType::dimension, GridImp, EntityImp>,
    //                               XT::Grid::extract_entity_t<EntityType, cd>>::value,
    //                  "Entity mismatch");
    CHECK_CRTP(this->as_imp().numDofs(entity));
    return this->as_imp().numDofs(entity);
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp().globalIndices(entity, ret));
  }

  Dune::DynamicVector<size_t> globalIndices(const EntityType& entity) const
  {
    Dune::DynamicVector<size_t> ret(numDofs(entity), 0);
    globalIndices(entity, ret);
    return ret;
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    CHECK_CRTP(this->as_imp(*this).mapToGlobal(entity, localIndex));
    return this->as_imp(*this).mapToGlobal(entity, localIndex);
  }
}; // class MapperInterface


class IsProductMapper
{
};

template <class Traits>
class ProductMapperInterface : public MapperInterface<Traits>, IsProductMapper
{
  typedef MapperInterface<Traits> BaseType;

public:
  using typename BaseType::EntityType;

  using BaseType::globalIndices;

  size_t size(const size_t factor_index) const
  {
    CHECK_CRTP(this->as_imp(*this).size(factor_index));
    return this->as_imp(*this).size(factor_index);
  }

  size_t numDofs(const size_t factor_index, const EntityType& entity) const
  {
    CHECK_CRTP(this->as_imp(*this).numDofs(factor_index, entity));
    return this->as_imp(*this).numDofs(factor_index, entity);
  }

  size_t maxNumDofs(const size_t factor_index) const
  {
    CHECK_CRTP(this->as_imp(*this).maxNumDofs(factor_index));
    return this->as_imp(*this).maxNumDofs(factor_index);
  }

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).globalIndices(factor_index, entity, ret));
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    CHECK_CRTP(this->as_imp(*this).mapToGlobal(factor_index, entity, local_index_in_factor));
    return this->as_imp(*this).mapToGlobal(factor_index, entity, local_index_in_factor);
  }

  size_t mapToLocal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    CHECK_CRTP(this->as_imp(*this).mapToLocal(factor_index, entity, local_index_in_factor));
    return this->as_imp(*this).mapToLocal(factor_index, entity, local_index_in_factor);
  }

  Dune::DynamicVector<size_t> globalIndices(const size_t factor_index, const EntityType& entity) const
  {
    Dune::DynamicVector<size_t> ret(numDofs(factor_index, entity), 0);
    globalIndices(factor_index, entity, ret);
    return ret;
  }
}; // class ProductMapperInterface


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_SPACES_MAPPER_INTERFACES_HH
