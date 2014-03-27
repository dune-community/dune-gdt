// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONAL_BASE_HH
#define DUNE_GDT_FUNCTIONAL_BASE_HH

#include <type_traits>

#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Functional {


// forward, to be used in the traits
template <class ImpTraits>
class AssemblableVolumeBase;


template <class ImpTraits>
class AssemblableVolumeBaseTraits
{
public:
  typedef typename ImpTraits::derived_type derived_type;
  typedef typename ImpTraits::GridViewType GridViewType;
  typedef typename ImpTraits::SpaceType SpaceType;
  typedef typename ImpTraits::VectorType VectorType;
  typedef typename ImpTraits::ScalarType ScalarType;

private:
  static_assert(std::is_base_of<Stuff::LA::VectorInterface<typename VectorType::Traits>, VectorType>::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");
}; // class AssemblableVolumeBaseTraits


template <class ImpTraits>
class AssemblableVolumeBase : public AssemblableFunctionalInterface<ImpTraits>,
                              public Functor::Codim0<typename ImpTraits::GridViewType>
{
  typedef AssemblableFunctionalInterface<ImpTraits> InterfaceType;
  typedef AssemblableVolumeBaseTraits<ImpTraits> Traits;
  typedef TmpStorageProvider::Vectors<typename ImpTraits::ScalarType> TmpStorageProviderType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::ScalarType ScalarType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;

private:
  typedef typename ImpTraits::LocalFunctionalType LocalFunctionalType;
  typedef LocalAssembler::Codim0Vector<LocalFunctionalType> LocalAssemblerType;

public:
  AssemblableVolumeBase(VectorType& vector, const SpaceType& space, const GridViewType& grid_view)
    : vector_(vector)
    , space_(space)
    , grid_view_(grid_view)
    , local_assembler_(nullptr)
    , tmp_storage_provider_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  AssemblableVolumeBase(VectorType& vector, const SpaceType& space)
    : vector_(vector)
    , space_(space)
    , grid_view_(*(space.grid_view()))
    , local_assembler_(nullptr)
    , tmp_storage_provider_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SpaceType& space() const
  {
    return space_;
  }

  VectorType& vector()
  {
    return vector_;
  }

  const VectorType& vector() const
  {
    return vector_;
  }

private:
  virtual const LocalFunctionalType& local_functional() const = 0;

public:
  virtual void prepare() DS_OVERRIDE
  {
    if (!assembled_ && !prepared_) {
      local_assembler_      = std::unique_ptr<LocalAssemblerType>(new LocalAssemblerType(local_functional()));
      tmp_storage_provider_ = std::unique_ptr<TmpStorageProviderType>(
          new TmpStorageProviderType(local_assembler_->numTmpObjectsRequired(), space_.mapper().maxNumDofs()));
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const EntityType& entity) DS_OVERRIDE
  {
    assert(prepared_);
    assert(local_assembler_);
    assert(tmp_storage_provider_);
    local_assembler_->assembleLocal(
        space_, entity, vector_, tmp_storage_provider_->vectors(), tmp_storage_provider_->indices());
  } // ... apply_local(...)

  void assemble()
  {
    if (!assembled_) {
      GridWalker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
      assembled_ = true;
    }
  } // ... assemble()

  template <class S>
  ScalarType apply(const Stuff::LA::VectorInterface<S>& source) const
  {
    typedef typename S::derived_type SourceType;
    assemble();
    return vector_.dot(static_cast<const SourceType&>(source));
  } // ... apply(...)

private:
  VectorType& vector_;
  const SpaceType& space_;
  const GridViewType& grid_view_;
  std::unique_ptr<LocalAssemblerType> local_assembler_;
  std::unique_ptr<TmpStorageProviderType> tmp_storage_provider_;
  bool prepared_;
  bool assembled_;
}; // class AssemblableVolumeBase


// forward, to be used in the traits
template <class ImpTraits>
class AssemblableFaceBase;


template <class ImpTraits>
class AssemblableFaceBaseTraits
{
public:
  typedef typename ImpTraits::derived_type derived_type;
  typedef typename ImpTraits::GridViewType GridViewType;
  typedef typename ImpTraits::SpaceType SpaceType;
  typedef typename ImpTraits::VectorType VectorType;
  typedef typename ImpTraits::ScalarType ScalarType;

private:
  static_assert(std::is_base_of<Stuff::LA::VectorInterface<typename VectorType::Traits>, VectorType>::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");
}; // class AssemblableFaceBaseTraits


template <class ImpTraits>
class AssemblableFaceBase : public AssemblableFunctionalInterface<ImpTraits>,
                            public Functor::Codim1<typename ImpTraits::GridViewType>
{
  typedef AssemblableFunctionalInterface<AssemblableFaceBaseTraits<ImpTraits>> InterfaceType;
  typedef AssemblableFaceBaseTraits<ImpTraits> Traits;
  typedef TmpStorageProvider::Vectors<typename ImpTraits::ScalarType> TmpStorageProviderType;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::ScalarType ScalarType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::Intersection IntersectionType;

private:
  typedef typename ImpTraits::LocalFunctionalType LocalFunctionalType;
  typedef LocalAssembler::Codim1Vector<LocalFunctionalType> LocalAssemblerType;

public:
  AssemblableFaceBase(VectorType& vector, const SpaceType& space, const GridViewType& grid_view)
    : vector_(vector)
    , space_(space)
    , grid_view_(grid_view)
    , local_assembler_(nullptr)
    , tmp_storage_provider_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  AssemblableFaceBase(VectorType& vector, const SpaceType& space)
    : vector_(vector)
    , space_(space)
    , grid_view_(*(space.grid_view()))
    , local_assembler_(nullptr)
    , tmp_storage_provider_(nullptr)
    , prepared_(false)
    , assembled_(false)
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SpaceType& space() const
  {
    return space_;
  }

  VectorType& vector()
  {
    return vector_;
  }

  const VectorType& vector() const
  {
    return vector_;
  }

private:
  virtual const LocalFunctionalType& local_functional() const = 0;

public:
  virtual void prepare() DS_OVERRIDE
  {
    if (!assembled_ && !prepared_) {
      local_assembler_      = std::unique_ptr<LocalAssemblerType>(new LocalAssemblerType(local_functional()));
      tmp_storage_provider_ = std::unique_ptr<TmpStorageProviderType>(
          new TmpStorageProviderType(local_assembler_->numTmpObjectsRequired(), space_.mapper().maxNumDofs()));
      prepared_ = true;
    }
  } // ... prepare()

  virtual void apply_local(const IntersectionType& intersection) DS_OVERRIDE
  {
    assert(prepared_);
    assert(local_assembler_);
    assert(tmp_storage_provider_);
    local_assembler_->assembleLocal(
        space_, intersection, vector_, tmp_storage_provider_->vectors(), tmp_storage_provider_->indices());
  } // ... apply_local(...)

  void assemble()
  {
    if (!assembled_) {
      GridWalker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
      assembled_ = true;
    }
  } // ... assemble()

  template <class S>
  ScalarType apply(const Stuff::LA::VectorInterface<S>& source) const
  {
    typedef typename S::derived_type SourceType;
    assemble();
    return vector_.dot(static_cast<const SourceType&>(source));
  } // ... apply(...)

private:
  VectorType& vector_;
  const SpaceType& space_;
  const GridViewType& grid_view_;
  std::unique_ptr<LocalAssemblerType> local_assembler_;
  std::unique_ptr<TmpStorageProviderType> tmp_storage_provider_;
  bool prepared_;
  bool assembled_;
}; // class AssemblableFaceBase


} // namespace Functional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONAL_BASE_HH
