// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#include <set>

#include <dune/geometry/type.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/xt/common/numeric_cast.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/local-functions.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/finite-elements/interfaces.hh>
#include <dune/gdt/spaces/mapper/interfaces.hh>

#ifndef DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH
#define DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH

namespace Dune {
namespace GDT {


template <class GL>
class ZeroOrderScalarDiscontinuousMapper;

template <class GL, class R = double, size_t r = 1, size_t rC = 1, class F = R>
class FixedOrderScalarDiscontinuousMapper;


namespace internal {


template <class GL>
class ZeroOrderScalarDiscontinuousMapperTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");

  template <int dim_>
  struct GeometryTypeLayout
  {
    bool contains(const GeometryType& gt) const
    {
      return gt.dim() == dim_;
    }
  };

public:
  using derived_type = ZeroOrderScalarDiscontinuousMapper<GL>;
  using BackendType = MultipleCodimMultipleGeomTypeMapper<GL, GeometryTypeLayout>;
  using EntityType = XT::Grid::extract_entity_t<GL>;
};


template <class GL, class R = double, size_t r = 1, size_t rC = 1, class F = R>
class FixedOrderScalarDiscontinuousMapperTraits
{
  static_assert(XT::Grid::is_layer<GL>::value, "");

  template <int dim_>
  struct GeometryTypeLayout
  {
    bool contains(const GeometryType& gt) const
    {
      return gt.dim() == dim_;
    }
  };

public:
  using derived_type = FixedOrderScalarDiscontinuousMapper<GL, R, r, rC, F>;
  using BackendType = MultipleCodimMultipleGeomTypeMapper<GL, GeometryTypeLayout>;
  using EntityType = XT::Grid::extract_entity_t<GL>;
};


} // namespace internal


template <class GL>
class ZeroOrderScalarDiscontinuousMapper
    : public MapperInterface<internal::ZeroOrderScalarDiscontinuousMapperTraits<GL>>
{
public:
  using Traits = internal::ZeroOrderScalarDiscontinuousMapperTraits<GL>;

private:
  using ThisType = ZeroOrderScalarDiscontinuousMapper<GL>;
  using BaseType = MapperInterface<Traits>;
  using D = typename GL::ctype;
  static const constexpr size_t d = GL::dimension;

public:
  using typename BaseType::EntityType;
  using typename BaseType::BackendType;

  ZeroOrderScalarDiscontinuousMapper(const GL& grid_layer)
    : mapper_(new BackendType(grid_layer))
  {
  }

  ZeroOrderScalarDiscontinuousMapper(const ThisType&) = default;
  ZeroOrderScalarDiscontinuousMapper(ThisType&&) = default;
  ZeroOrderScalarDiscontinuousMapper& operator=(const ThisType&) = delete;
  ZeroOrderScalarDiscontinuousMapper& operator=(ThisType&&) = delete;

  const BackendType& backend() const
  {
    return *mapper_;
  }

  size_t size() const
  {
    return mapper_->size();
  }

  size_t maxNumDofs() const
  {
    return 1;
  }

  size_t numDofs(const EntityType& /*entity*/) const
  {
    return 1;
  }

  template <int cd, class GridImp, template <int, int, class> class EntityImp>
  typename std::enable_if<cd != EntityType::codimension, size_t>::type
  numDofs(const Entity<cd, EntityType::dimension, GridImp, EntityImp>& /*entity*/) const
  {
    return 0;
  }

  using BaseType::globalIndices;

  void globalIndices(const EntityType& entity, DynamicVector<size_t>& ret) const
  {
    ret[0] = mapper_->index(entity);
  }

  size_t mapToGlobal(const EntityType& entity, const size_t local_index) const
  {
    if (local_index >= 1)
      DUNE_THROW(Exception, "numDofs(entity) = " << numDofs(entity) << "\n   local_index = " << local_index);
    return mapper_->index(entity);
  }

private:
  const std::shared_ptr<BackendType> mapper_;
}; // class ZeroOrderScalarDiscontinuousMapper


template <class GL, class R, size_t r, size_t rC, class F>
class FixedOrderScalarDiscontinuousMapper
    : public MapperInterface<internal::FixedOrderScalarDiscontinuousMapperTraits<GL, R, r, rC, F>>
{
public:
  using Traits = internal::FixedOrderScalarDiscontinuousMapperTraits<GL, R, r, rC, F>;

private:
  using ThisType = FixedOrderScalarDiscontinuousMapper<GL, R, r, rC, F>;
  using BaseType = MapperInterface<Traits>;
  using D = typename GL::ctype;
  static const constexpr size_t d = GL::dimension;

public:
  using FiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC, F>;
  using typename BaseType::EntityType;
  using typename BaseType::BackendType;

  FixedOrderScalarDiscontinuousMapper(
      const GL& grid_layer,
      const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>>& finite_elements)
    : finite_elements_(finite_elements)
    , mapper_(new BackendType(grid_layer))
    , offset_(new std::vector<size_t>(mapper_->size()))
    , max_num_dofs_(0)
    , size_(0)
  {
    for (auto&& element : elements(grid_layer)) {
      (*offset_)[mapper_->subIndex(element, 0, 0)] = size_;
      const auto finite_element_search_result = finite_elements_->find(element.geometry().type());
      if (finite_element_search_result == finite_elements_->end())
        DUNE_THROW(mapper_error,
                   "Missing finite element for the required geometry type " << element.geometry().type() << "!");
      const auto& finite_element = *finite_element_search_result->second;
      size_ += finite_element.size();
      max_num_dofs_ = std::max(max_num_dofs_, XT::Common::numeric_cast<size_t>(finite_element.size()));
    }
  } // ... FixedOrderScalarDiscontinuousMapper(...)

  FixedOrderScalarDiscontinuousMapper(const ThisType&) = default;
  FixedOrderScalarDiscontinuousMapper(ThisType&&) = default;
  FixedOrderScalarDiscontinuousMapper& operator=(const ThisType&) = delete;
  FixedOrderScalarDiscontinuousMapper& operator=(ThisType&&) = delete;

  const BackendType& backend() const
  {
    return *mapper_;
  }

  size_t size() const
  {
    return size_;
  }

  size_t maxNumDofs() const
  {
    return max_num_dofs_;
  }

  size_t numDofs(const EntityType& entity) const
  {
    const auto finite_element_search_result = finite_elements_->find(entity.geometry().type());
    if (finite_element_search_result == finite_elements_->end())
      DUNE_THROW(XT::Common::Exceptions::internal_error,
                 "This must not happen after the checks in the ctor, the grid layer did not report all geometry types!"
                     << "\n   entity.geometry().type() = "
                     << entity.geometry().type());
    const auto& finite_element = *finite_element_search_result->second;
    return finite_element.size();
  } // ... numDofs(...)

  template <int cd, class GridImp, template <int, int, class> class EntityImp>
  typename std::enable_if<cd != EntityType::codimension, size_t>::type
  numDofs(const Entity<cd, EntityType::dimension, GridImp, EntityImp>& /*entity*/) const
  {
    return 0;
  }

  using BaseType::globalIndices;

  void globalIndices(const EntityType& entity, DynamicVector<size_t>& ret) const
  {
    const size_t offset = (*offset_)[mapper_->subIndex(entity, 0, 0)];
    const auto local_size = numDofs(entity);
    if (ret.size() < local_size)
      ret.resize(local_size, 0);
    for (size_t ii = 0; ii < local_size; ++ii)
      ret[ii] = offset + ii;
  } // ... globalIndices(...)

  size_t mapToGlobal(const EntityType& entity, const size_t local_index) const
  {
    if (local_index >= numDofs(entity))
      DUNE_THROW(Exception, "numDofs(entity) = " << numDofs(entity) << "\n   local_index = " << local_index);
    const size_t offset = (*offset_)[mapper_->subIndex(entity, 0, 0)];
    return offset + local_index;
  }

private:
  const std::shared_ptr<std::map<GeometryType, std::shared_ptr<FiniteElementType>>> finite_elements_;
  const std::shared_ptr<BackendType> mapper_;
  std::shared_ptr<std::vector<size_t>> offset_;
  size_t max_num_dofs_;
  size_t size_;
}; // class FixedOrderScalarDiscontinuousMapper


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_DISCONTINUOUS_HH
