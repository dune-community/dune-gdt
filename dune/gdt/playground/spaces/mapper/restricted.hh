// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_SPACES_MAPPER_RESTRICTED_HH
#define DUNE_GDT_SPACES_MAPPER_RESTRICTED_HH

#include <dune/grid/common/rangegenerators.hh>

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


// forward
template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedMapper;


namespace internal {


template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedMapperTraits
{
  static_assert(XT::Grid::is_layer<RestrictionGridLayer>::value, "");
  static_assert(is_space<UnrestrictedSpace>::value, "");

public:
  typedef RestrictedMapper<UnrestrictedSpace, RestrictionGridLayer> derived_type;
  typedef typename UnrestrictedSpace::MapperType BackendType;
  typedef XT::Grid::extract_entity_t<RestrictionGridLayer> EntityType;
};


} // namespace internal


template <class UnrestrictedSpace, class RestrictionGridLayer>
class RestrictedMapper
    : public MapperInterface<internal::RestrictedMapperTraits<UnrestrictedSpace, RestrictionGridLayer>>
{
  typedef MapperInterface<internal::RestrictedMapperTraits<UnrestrictedSpace, RestrictionGridLayer>> BaseType;
  typedef RestrictedMapper<UnrestrictedSpace, RestrictionGridLayer> ThisType;

public:
  using typename BaseType::BackendType;
  using typename BaseType::EntityType;

  RestrictedMapper(const UnrestrictedSpace& unrestricted_space, RestrictionGridLayer restriction_grid_layer)
    : unrestricted_space_(unrestricted_space)
    , grid_layer_(restriction_grid_layer)
    , max_num_dofs_(0)
    , map_to_restricted_()
    , map_to_unrestricted_()
  {
    std::set<size_t> restricted_indices;
    DynamicVector<size_t> unrestricted_indices(backend().maxNumDofs(), 0);
    for (auto&& entity : elements(grid_layer_)) {
      max_num_dofs_ = std::max(max_num_dofs_, backend().numDofs(entity));
      backend().globalIndices(entity, unrestricted_indices);
      for (size_t ii = 0; ii < backend().numDofs(entity); ++ii)
        restricted_indices.insert(unrestricted_indices[ii]);
    }
    map_to_unrestricted_ = std::vector<size_t>(restricted_indices.size());
    size_t counter = 0;
    for (const auto& index : restricted_indices) {
      map_to_restricted_[index] = counter;
      map_to_unrestricted_[counter] = index;
      ++counter;
    }
    if (map_to_restricted_.size() != map_to_unrestricted_.size())
      DUNE_THROW(InvalidStateException, "This should not happen!");
  }

  RestrictedMapper(const ThisType& other) = default;
  RestrictedMapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  template <class U, class V>
  void restrict(const XT::LA::VectorInterface<U>& unrestricted_vector,
                XT::LA::VectorInterface<V>& restricted_vector) const
  {
    if (unrestricted_vector.size() != backend().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "unrestricted_vector.size() = " << unrestricted_vector.size()
                                                 << "\n   unrestricted_space.mapper().size() = "
                                                 << backend().size());
    if (restricted_vector.size() != size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "restricted_vector.size() = " << restricted_vector.size() << "\n   size() = " << size());
    // the actual work
    for (size_t restricted_DoF = 0; restricted_DoF < map_to_unrestricted_.size(); ++restricted_DoF)
      restricted_vector[restricted_DoF] = unrestricted_vector[map_to_unrestricted_[restricted_DoF]];
  } // ... restrict(...)

  template <class V>
  typename V::derived_type restrict(const XT::LA::VectorInterface<V>& unrestricted_vector) const
  {
    typename V::derived_type restricted_vector(size(), 0.);
    restrict(unrestricted_vector, restricted_vector);
    return restricted_vector;
  }

  template <class U, class V>
  void extend(const XT::LA::VectorInterface<U>& restricted_vector,
              XT::LA::VectorInterface<V>& unrestricted_vector) const
  {
    if (restricted_vector.size() != size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "restricted_vector.size() = " << restricted_vector.size() << "\n   size() = " << size());
    if (unrestricted_vector.size() != backend().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "unrestricted_vector.size() = " << unrestricted_vector.size()
                                                 << "\n   unrestricted_space.mapper().size() = "
                                                 << backend().size());
    // the actual work
    unrestricted_vector *= 0.;
    for (size_t restricted_DoF = 0; restricted_DoF < map_to_unrestricted_.size(); ++restricted_DoF)
      unrestricted_vector[map_to_unrestricted_[restricted_DoF]] = restricted_vector[restricted_DoF];
  } // ... extend(...)

  template <class V>
  typename V::derived_type extend(const XT::LA::VectorInterface<V>& restricted_vector) const
  {
    typename V::derived_type unrestricted_vector(backend().size(), 0.);
    extend(restricted_vector, unrestricted_vector);
    return unrestricted_vector;
  }

  const BackendType& backend() const
  {
    return unrestricted_space_.mapper();
  }

  size_t size() const
  {
    return map_to_restricted_.size();
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  size_t numDofs(const EntityType& entity) const
  {
    const auto error_message = check_entity(entity);
    if (error_message.size() > 0)
      DUNE_THROW(restricted_space_error, error_message);
    return backend().numDofs(entity);
  }

  void globalIndices(const EntityType& entity, DynamicVector<size_t>& ret) const
  {
    const auto error_message = check_entity(entity);
    if (error_message.size() > 0)
      DUNE_THROW(restricted_space_error, error_message);
    backend().globalIndices(entity, ret);
    for (size_t ii = 0; ii < numDofs(entity); ++ii)
      ret[ii] = map_to_restricted_.at(ret[ii]);
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    const auto error_message = check_entity(entity);
    if (error_message.size() > 0)
      DUNE_THROW(restricted_space_error, error_message);
    return map_to_restricted_.at(backend().mapToGlobal(entity, localIndex));
  }

private:
  std::string check_entity(const EntityType& entity) const
  {
    std::stringstream error_message;
    if (!grid_layer_.indexSet().contains(entity)) {
      if (unrestricted_space_.grid_layer().indexSet().contains(entity))
        error_message << "Entity not contained in restriction grid layer, but contained in the unrestricted grid layer "
                         "with index "
                      << unrestricted_space_.grid_layer().indexSet().index(entity) << "!";
      else
        error_message << "Entity neither contained in restriction grid layer nor in the unrestricted grid layer!";
    }
    return error_message.str();
  } // ... check_entity(...)

  const UnrestrictedSpace unrestricted_space_;
  const RestrictionGridLayer grid_layer_;
  size_t max_num_dofs_;
  std::map<size_t, size_t> map_to_restricted_;
  std::vector<size_t> map_to_unrestricted_;
}; // class RestrictedMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_RESTRICTED_HH
