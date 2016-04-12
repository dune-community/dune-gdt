// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

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
namespace Mapper {

#if HAVE_DUNE_PDELAB


// forwards
template <class PdelabSpaceImp, size_t rangeDim>
class ContinuousPdelabWrapper;

template <class PdelabSpaceImp>
class DiscontinuousPdelabWrapper;


namespace internal {


template <class PdelabSpaceImp, size_t rangeDim>
class ContinuousPdelabWrapperTraits
{
public:
  typedef ContinuousPdelabWrapper<PdelabSpaceImp, rangeDim> derived_type;
  typedef PdelabSpaceImp BackendType;
  typedef typename BackendType::Element EntityType;
};

template <class PdelabSpaceImp>
class DiscontinuousPdelabWrapperTraits
{
public:
  typedef DiscontinuousPdelabWrapper<PdelabSpaceImp> derived_type;
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
class ContinuousPdelabWrapper
    : public DefaultProductMapperFromTuple<
          typename PdelabSpaceImp::Traits::GridViewType,
          typename DSC::make_identical_tuple<ContinuousPdelabWrapper<PdelabSpaceImp, 1>, rangeDim>::type>::type
{
  typedef ContinuousPdelabWrapper<PdelabSpaceImp, 1> ScalarValuedMapperType;
  typedef typename DefaultProductMapperFromTuple<
      typename PdelabSpaceImp::Traits::GridViewType,
      typename Dune::Stuff::Common::make_identical_tuple<ScalarValuedMapperType, rangeDim>::type>::type BaseType;

public:
  typedef typename internal::ContinuousPdelabWrapperTraits<PdelabSpaceImp, rangeDim>::BackendType BackendType;
  ContinuousPdelabWrapper(const BackendType& pdelab_space)
    : BaseType(pdelab_space.gridView(),
               DSC::make_identical_tuple<ScalarValuedMapperType, rangeDim>::create(pdelab_space))
  {
  }
}; // class ContinuousPdelabWrapper


template <class PdelabSpaceImp>
class ContinuousPdelabWrapper<PdelabSpaceImp, 1>
    : public internal::PdelabWrapperBase<internal::ContinuousPdelabWrapperTraits<PdelabSpaceImp, 1>>
{
public:
  typedef typename internal::ContinuousPdelabWrapperTraits<PdelabSpaceImp, 1> Traits;
  typedef typename Traits::EntityType EntityType;

  template <class... Args>
  ContinuousPdelabWrapper(Args&&... args)
    : internal::PdelabWrapperBase<Traits>(std::forward<Args>(args)...)
  {
  }

protected:
  virtual size_t mapAfterBound(const EntityType& /*entity*/, const size_t& localIndex) const override
  {
    return this->lfs_.dofIndex(localIndex).entityIndex()[1];
  }
}; // class ContinuousPdelabWrapper


template <class PdelabSpaceImp>
class DiscontinuousPdelabWrapper
    : public internal::PdelabWrapperBase<internal::DiscontinuousPdelabWrapperTraits<PdelabSpaceImp>>
{
public:
  typedef typename internal::DiscontinuousPdelabWrapperTraits<PdelabSpaceImp> Traits;
  typedef typename Traits::EntityType EntityType;

  template <class... Args>
  DiscontinuousPdelabWrapper(Args&&... args)
    : internal::PdelabWrapperBase<Traits>(std::forward<Args>(args)...)
  {
  }

protected:
  virtual size_t mapAfterBound(const EntityType& entity, const size_t& localIndex) const override
  {
    return this->lfs_.dofIndex(localIndex).entityIndex()[1] * this->numDofs(entity) + localIndex;
  }
}; // class DiscontinuousPdelabWrapper


#else // HAVE_DUNE_PDELAB


template <class PdelabSpaceImp>
class ContinuousPdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};

template <class PdelabSpaceImp>
class DiscontinuousPdelabWrapper
{
  static_assert(Dune::AlwaysFalse<PdelabSpaceImp>::value, "You are missing dune-pdelab!");
};


#endif // HAVE_DUNE_PDELAB

} // namespace Mapper
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_MAPPER_DUNE_PDELAB_WRAPPER_HH
