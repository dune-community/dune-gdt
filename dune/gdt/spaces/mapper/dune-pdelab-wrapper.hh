// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014, 2016)

#ifndef DUNE_GDT_SPACES_MAPPER_DUNE_PDELAB_WRAPPER_HH
#define DUNE_GDT_SPACES_MAPPER_DUNE_PDELAB_WRAPPER_HH

#include <unordered_map>

#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>

#if HAVE_DUNE_PDELAB
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#endif

#include <dune/stuff/common/parallel/threadstorage.hh>
#include <dune/stuff/common/type_utils.hh>
#include <dune/stuff/common/tuple.hh>

#include "interfaces.hh"
#include "product.hh"

namespace Dune {
namespace GDT {

#if HAVE_DUNE_PDELAB


// forwards
template <class PdelabSpaceImp, size_t rangeDim>
class DunePdelabCgMapperWrapper;

template <class PdelabSpaceImp>
class DunePdelabDgMapperWrapper;


namespace internal {


template <class PdelabSpaceImp, size_t rangeDim>
class DunePdelabCgMapperWrapperTraits
{
public:
  typedef DunePdelabCgMapperWrapper<PdelabSpaceImp, rangeDim> derived_type;
  typedef PdelabSpaceImp BackendType;
  typedef typename BackendType::Element EntityType;
};

template <class PdelabSpaceImp>
class DunePdelabDgMapperWrapperTraits
{
public:
  typedef DunePdelabDgMapperWrapper<PdelabSpaceImp> derived_type;
  typedef PdelabSpaceImp BackendType;
  typedef typename BackendType::Element EntityType;
};


template <class ImpTraits>
class PdelabWrapperBase : public MapperInterface<ImpTraits>
{
  typedef MapperInterface<ImpTraits> InterfaceType;

public:
  typedef typename InterfaceType::EntityType EntityType;
  typedef typename InterfaceType::BackendType BackendType;

private:
  typedef PDELab::LocalFunctionSpace<BackendType, PDELab::TrialSpaceTag> PdeLabLFSType;

public:
  typedef typename PdeLabLFSType::Traits::DOFIndex MultiIndexType;

public:
  explicit PdelabWrapperBase(const BackendType& pdelab_space)
    : backend_(pdelab_space)
    , lfs_(backend_)
  {
    const auto& grid_view = backend_.gridView();
    const auto it_end     = grid_view.template end<0>();
    std::size_t count = 0;
    for (auto it = grid_view.template begin<0>(); it != it_end; ++it) {
      const auto& entity = *it;
      lfs_.bind(entity);
      for (size_t ii = 0; ii < lfs_.size(); ++ii)
        if (index_map_.find(lfs_.dofIndex(ii)) == index_map_.end())
          index_map_.insert(std::make_pair(lfs_.dofIndex(ii), count++));
    }
  }

  virtual ~PdelabWrapperBase()
  {
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  size_t numDofs(const EntityType& entity) const
  {
    lfs_.bind(entity);
    return lfs_.size();
  }

  size_t maxNumDofs() const
  {
    return backend_.maxLocalSize();
  }

  void globalIndices(const EntityType& entity, Dune::DynamicVector<size_t>& ret) const
  {
    lfs_.bind(entity);
    // some checks
    const size_t numLocalDofs = numDofs(entity);
    if (ret.size() < numLocalDofs)
      ret.resize(numLocalDofs);
    // compute
    for (size_t ii = 0; ii < numLocalDofs; ++ii)
      ret[ii] = mapToGlobal(entity, ii);
  } // ... globalIndices(...)

  using InterfaceType::globalIndices;

  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    lfs_.bind(entity);
    assert(localIndex < lfs_.size());
    return index_map_[lfs_.dofIndex(localIndex)];
  } // ... mapToGlobal(...)

protected:
  virtual size_t mapAfterBound(const EntityType& entity, const size_t& localIndex) const = 0;

  const BackendType& backend_;
  mutable PdeLabLFSType lfs_;
  mutable std::unordered_map<MultiIndexType, std::size_t> index_map_;
}; // class PdelabWrapperBase


} // namespace internal


template <class PdelabSpaceImp, size_t rangeDim = 1>
class DunePdelabCgMapperWrapper
    : public DefaultProductMapperFromTuple<
          typename PdelabSpaceImp::Traits::GridViewType,
          typename DSC::make_identical_tuple<DunePdelabCgMapperWrapper<PdelabSpaceImp, 1>, rangeDim>::type>::type
{
  typedef DunePdelabCgMapperWrapper<PdelabSpaceImp, 1> ScalarValuedMapperType;
  typedef typename DefaultProductMapperFromTuple<
      typename PdelabSpaceImp::Traits::GridViewType,
      typename Dune::Stuff::Common::make_identical_tuple<ScalarValuedMapperType, rangeDim>::type>::type BaseType;

public:
  typedef typename internal::DunePdelabCgMapperWrapperTraits<PdelabSpaceImp, rangeDim>::BackendType BackendType;
  DunePdelabCgMapperWrapper(const BackendType& pdelab_space)
    : BaseType(pdelab_space.gridView(),
               DSC::make_identical_tuple<ScalarValuedMapperType, rangeDim>::create(pdelab_space))
  {
  }
}; // class DunePdelabCgMapperWrapper


template <class PdelabSpaceImp>
class DunePdelabCgMapperWrapper<PdelabSpaceImp, 1>
    : public internal::PdelabWrapperBase<internal::DunePdelabCgMapperWrapperTraits<PdelabSpaceImp, 1>>
{
public:
  typedef typename internal::DunePdelabCgMapperWrapperTraits<PdelabSpaceImp, 1> Traits;
  typedef typename Traits::EntityType EntityType;

  template <class... Args>
  DunePdelabCgMapperWrapper(Args&&... args)
    : internal::PdelabWrapperBase<Traits>(std::forward<Args>(args)...)
  {
  }

protected:
  virtual size_t mapAfterBound(const EntityType& /*entity*/, const size_t& localIndex) const override
  {
    return this->lfs_.dofIndex(localIndex).entityIndex()[1];
  }
}; // class DunePdelabCgMapperWrapper


template <class PdelabSpaceImp>
class DunePdelabDgMapperWrapper
    : public internal::PdelabWrapperBase<internal::DunePdelabDgMapperWrapperTraits<PdelabSpaceImp>>
{
public:
  typedef typename internal::DunePdelabDgMapperWrapperTraits<PdelabSpaceImp> Traits;
  typedef typename Traits::EntityType EntityType;

  template <class... Args>
  DunePdelabDgMapperWrapper(Args&&... args)
    : internal::PdelabWrapperBase<Traits>(std::forward<Args>(args)...)
  {
  }

protected:
  virtual size_t mapAfterBound(const EntityType& entity, const size_t& localIndex) const override
  {
    return this->lfs_.dofIndex(localIndex).entityIndex()[1] * this->numDofs(entity) + localIndex;
  }
}; // class DunePdelabDgMapperWrapper


#else // HAVE_DUNE_PDELAB


template <class PdelabSpaceImp>
class DunePdelabCgMapperWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};

template <class PdelabSpaceImp>
class DunePdelabDgMapperWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_DUNE_PDELAB_WRAPPER_HH
