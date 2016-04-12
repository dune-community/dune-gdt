// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH
#define DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH

#include <dune/common/dynvector.hh>

#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/spaces/mapper/interfaces.hh>
#include <dune/gdt/spaces/mapper/dune-pdelab-wrapper.hh>

namespace Dune {
namespace GDT {
namespace Mapper {


// forward
template <class PdelabSpaceImp, size_t rangeDim = 1, size_t rangeDimCols = 1>
class ProductDG
{
  static_assert(AlwaysFalse<PdelabSpaceImp>::value, "Not available for these dimensions!");
};


namespace internal {


template <class PdelabSpaceImp, size_t rangeDim, size_t rangeDimCols>
class ProductDGTraits
{
  static_assert(rangeDim >= 1, "Really?");
  static_assert(rangeDimCols >= 1, "Really?");

public:
  typedef ProductDG<PdelabSpaceImp, rangeDim, rangeDimCols> derived_type;
  typedef PdelabSpaceImp BackendType;
  typedef typename BackendType::Element EntityType;
};


} // namespace internal


template <class PdelabSpaceImp, size_t rangeDim>
class ProductDG<PdelabSpaceImp, rangeDim, 1>
    : public ProductMapperInterface<internal::ProductDGTraits<PdelabSpaceImp, rangeDim, 1>>
{
  typedef ProductMapperInterface<internal::ProductDGTraits<PdelabSpaceImp, rangeDim, 1>> InterfaceType;
  static const size_t dimRange = rangeDim;

public:
  typedef internal::ProductDGTraits<PdelabSpaceImp, rangeDim, 1> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef DiscontinuousPdelabWrapper<BackendType> FactorMapperType;
  typedef typename Traits::EntityType EntityType;

  ProductDG(const BackendType& pdelab_space)
    : backend_(pdelab_space)
    , factor_mapper_(pdelab_space)
  {
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return factor_mapper_.size() * dimRange;
  }

  size_t numDofs(const EntityType& entity) const
  {
    return dimRange * factor_mapper_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return dimRange * factor_mapper_.maxNumDofs();
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    if (ret.size() < numDofs(entity))
      ret.resize(numDofs(entity));
    const auto factor_num_dofs       = factor_mapper_.numDofs(entity);
    const auto factor_mapper_size    = factor_mapper_.size();
    const auto factor_global_indices = factor_mapper_.globalIndices(entity);
    for (size_t ii = 0; ii < dimRange; ++ii) {
      const auto factor_num_dofs_times_ii    = factor_num_dofs * ii;
      const auto factor_mapper_size_times_ii = factor_mapper_size * ii;
      for (size_t jj = 0; jj < factor_num_dofs; ++jj) {
        ret[factor_num_dofs_times_ii + jj] = factor_global_indices[jj] + factor_mapper_size_times_ii;
      }
    }
  } // ... globalIndices(...)

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < numDofs(entity));
    size_t factor_index               = 0;
    const auto factor_mapper_num_dofs = factor_mapper_.numDofs(entity);
    while (localIndex >= factor_mapper_num_dofs) {
      localIndex -= factor_mapper_num_dofs;
      ++factor_index;
    }
    return factor_mapper_.mapToGlobal(entity, localIndex) + factor_index * factor_mapper_.size();
  }

  using InterfaceType::globalIndices;

  size_t numDofs(const size_t /*factor_index*/, const EntityType& entity) const
  {
    return factor_mapper_.numDofs(entity);
  }

  size_t maxNumDofs(const size_t /*factor_index*/) const
  {
    return factor_mapper_.maxNumDofs();
  }

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    assert(factor_index < dimRange);
    const auto factor_mapper_num_dofs = factor_mapper_.numDofs(entity);
    if (ret.size() < factor_mapper_num_dofs)
      ret.resize(factor_mapper_num_dofs);
    const auto factor_mapper_size_times_factor_index = factor_mapper_.size() * factor_index;
    const auto factor_mapper_global_indices = factor_mapper_.globalIndices(entity);
    for (size_t jj = 0; jj < factor_mapper_num_dofs; ++jj)
      ret[jj] = factor_mapper_global_indices[jj] + factor_mapper_size_times_factor_index;
  } // ... globalIndices(...)

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    assert(local_index_in_factor < factor_mapper_.numDofs(entity));
    assert(factor_index < dimRange);
    return factor_mapper_.mapToGlobal(entity, local_index_in_factor) + factor_index * factor_mapper_.size();
  }

  size_t mapToLocal(const size_t factor_index, const EntityType& entity, const size_t& local_index_in_factor) const
  {
    assert(local_index_in_factor < factor_mapper_.numDofs(entity));
    assert(factor_index < dimRange);
    return factor_mapper_.numDofs(entity) * factor_index + local_index_in_factor;
  }

private:
  const BackendType& backend_;
  const FactorMapperType factor_mapper_;
}; // class ProductDG< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH
