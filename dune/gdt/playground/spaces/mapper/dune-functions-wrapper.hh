// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_PLAYGROUND_SPACES_MAPPER_DUNE_FUNCTIONS_WRAPPER_HH
#define DUNE_GDT_PLAYGROUND_SPACES_MAPPER_DUNE_FUNCTIONS_WRAPPER_HH

#include <dune/common/typetraits.hh>

#if HAVE_DUNE_FUNCTIONS
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#endif

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/spaces/mapper/interfaces.hh>

namespace Dune {
namespace GDT {

#if HAVE_DUNE_FUNCTIONS


// forward, to be used in the traits and to allow for specialization
template <class GL, int p, class R, size_t r, size_t rC>
class DuneFunctionsMapperWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "Untested for these dimensions!");
};


namespace internal {


template <class GL, int p, class R, size_t r, size_t rC>
class DuneFunctionsMapperWrapperTraits
{
  static_assert(XT::Grid::is_view<GL>::value, "We Probably need to use TemporaryGridView from dune-xt-grid!");

public:
  typedef DuneFunctionsMapperWrapper<GL, p, R, r, rC> derived_type;
  typedef Functions::LagrangeDGBasis<GL, p> BackendType;
  typedef typename XT::Grid::extract_entity<GL>::type EntityType;
};


} // namespace internal


template <class GL, int p, class R>
class DuneFunctionsMapperWrapper<GL, p, R, 1, 1>
    : public MapperInterface<internal::DuneFunctionsMapperWrapperTraits<GL, p, R, 1, 1>>
{
  typedef DuneFunctionsMapperWrapper<GL, p, R, 1, 1> ThisType;
  typedef MapperInterface<internal::DuneFunctionsMapperWrapperTraits<GL, p, R, 1, 1>> BaseType;

public:
  typedef internal::DuneFunctionsMapperWrapperTraits<GL, p, R, 1, 1> Traits;

  using typename BaseType::BackendType;
  using typename BaseType::EntityType;

  DuneFunctionsMapperWrapper(std::shared_ptr<const BackendType> bcknd)
    : backend_(bcknd)
    , max_num_dofs_(0)
  {
    auto local_view = backend_->localView();
    auto local_index_set = backend_->localIndexSet();
    for (auto&& entity : elements(backend_->gridView())) {
      local_view.bind(entity);
      local_index_set.bind(local_view);
      max_num_dofs_ = std::max(max_num_dofs_, local_index_set.size());
    }
  }

  DuneFunctionsMapperWrapper(const ThisType& other) = default;
  DuneFunctionsMapperWrapper(ThisType&& source) = default;

  ThisType& operator=(const ThisType& other) = delete;
  ThisType& operator=(ThisType&& source) = delete;

  const BackendType& backend() const
  {
    return *backend_;
  }

  size_t size() const
  {
    return backend_->size();
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  size_t numDofs(const EntityType& entity) const
  {
    auto local_view = backend_->localView();
    auto local_index_set = backend_->localIndexSet();
    local_view.bind(entity);
    local_index_set.bind(local_view);
    return local_index_set.size();
  }

  using BaseType::globalIndices;

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    auto local_view = backend_->localView();
    auto local_index_set = backend_->localIndexSet();
    local_view.bind(entity);
    local_index_set.bind(local_view);
    assert(ret.size() >= local_index_set.size());
    for (size_t ii = 0; ii < local_index_set.size(); ++ii)
      ret[ii] = local_index_set.index(ii)[0];
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    auto local_view = backend_->localView();
    auto local_index_set = backend_->localIndexSet();
    local_view.bind(entity);
    local_index_set.bind(local_view);
    assert(localIndex < local_index_set.size());
    return local_index_set.index(localIndex)[0];
  }

private:
  const std::shared_ptr<const BackendType> backend_;
  size_t max_num_dofs_;
}; // class DuneFunctionsMapperWrapper


#else // HAVE_DUNE_FUNCTIONS


template <class GL, int p, class R, size_t r, size_t rC = 1>
class DuneFunctionsMapperWrapper
{
  static_assert(Dune::AlwaysFalse<GL>::value, "You are missing dune-functions!");
};


#endif // HAVE_DUNE_FUNCTIONS

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_SPACES_MAPPER_DUNE_FUNCTIONS_WRAPPER_HH
