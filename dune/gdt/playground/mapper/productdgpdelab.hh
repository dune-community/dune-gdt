// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH
#define DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH

#include <dune/common/dynvector.hh>

#include <dune/stuff/common/debug.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/mapper/interface.hh>
#include <dune/gdt/mapper/pdelab.hh>

namespace Dune {
namespace GDT {
namespace Mapper {


// forward
template< class PdelabSpaceImp, size_t rangeDim = 1, size_t rangeDimCols = 1 >
class ProductDG
{
  static_assert(AlwaysFalse< PdelabSpaceImp >::value, "Not available for these dimensions!");
};


namespace internal {


template< class PdelabSpaceImp, size_t rangeDim, size_t rangeDimCols >
class ProductDGTraits
{
  static_assert(rangeDim >= 1, "Really?");
  static_assert(rangeDimCols >= 1, "Really?");
public:
  typedef ProductDG< PdelabSpaceImp, rangeDim, rangeDimCols> derived_type;
  typedef PdelabSpaceImp                                     BackendType;
  typedef typename BackendType::Element                      EntityType;
};


} // namespace internal


template< class PdelabSpaceImp, size_t rangeDim >
class ProductDG< PdelabSpaceImp, rangeDim, 1 >
  : public MapperInterface< internal::ProductDGTraits< PdelabSpaceImp, rangeDim, 1 > >
{
  typedef MapperInterface< internal::ProductDGTraits< PdelabSpaceImp, rangeDim, 1 > > InterfaceType;
  static const size_t dimRange = rangeDim;
public:
  typedef internal::ProductDGTraits< PdelabSpaceImp, rangeDim, 1 > Traits;
  typedef typename Traits::BackendType                             BackendType;
  typedef DiscontinuousPdelabWrapper< BackendType >                FactorMapperType;
  typedef typename Traits::EntityType                              EntityType;

  ProductDG(const BackendType& pdelab_space)
    : backend_(pdelab_space)
    , factor_mapper_(pdelab_space)
  {}

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return factor_mapper_.size()*dimRange;
  }

  size_t numDofs(const EntityType& entity) const
  {
    return dimRange*factor_mapper_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return dimRange*factor_mapper_.maxNumDofs();
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    if (ret.size() < numDofs(entity))
      ret.resize(numDofs(entity));
    for (size_t ii = 0; ii < dimRange; ++ii) {
      for (size_t jj = 0; jj < factor_mapper_.numDofs(entity); ++jj) {
        ret[ii*factor_mapper_.numDofs(entity)+jj] = factor_mapper_.globalIndices(entity)[jj] + ii*factor_mapper_.size();
      }
    }
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  void globalIndices(const size_t factor_index, const EntityType& entity, Dune::DynamicVector< size_t >& ret) const
  {
    assert(factor_index < dimRange);
    if (ret.size() < factor_mapper_.numDofs(entity))
      ret.resize(factor_mapper_.numDofs(entity));
    for (size_t jj = 0; jj < factor_mapper_.numDofs(entity); ++jj)
      ret[jj] = factor_mapper_.globalIndices(entity)[jj] + factor_index*factor_mapper_.size();
  } // ... globalIndices(...)

  Dune::DynamicVector< size_t > globalIndices(const size_t factor_index, const EntityType& entity) const
  {
    Dune::DynamicVector< size_t > ret(factor_mapper_.numDofs(entity), 0);
    globalIndices(factor_index, entity, ret);
    return ret;
  }

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < numDofs(entity));
    size_t factor_index = 0;
    while (localIndex >= factor_mapper_.numDofs(entity)) {
      localIndex -= factor_mapper_.numDofs(entity);
      ++factor_index;
    }
    return factor_mapper_.globalIndices(entity)[localIndex] + factor_index*factor_mapper_.size();
  }

  size_t mapToGlobal(const size_t factor_index, const EntityType& entity, const size_t& localIndex) const
  {
    assert(localIndex < factor_mapper_.numDofs(entity));
    assert(factor_index < dimRange);
    return factor_mapper_.globalIndices(entity)[localIndex] + factor_index*factor_mapper_.size();
  }

private:
  const BackendType& backend_;
  const FactorMapperType factor_mapper_;
}; // class ProductFiniteVolume< ..., rangeDim, 1 >


} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_MAPPER_PRODUCTDGPDELAB_HH
