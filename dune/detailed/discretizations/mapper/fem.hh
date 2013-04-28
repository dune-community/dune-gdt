#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH

#include <dune/common/dynvector.hh>

#include <dune/fem/space/mapper/dofmapper.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Mapper {


// forward, to be used in the traits and to allow for specialization
template< class SpaceTraits >
class FemWrapper;


template< class SpaceTraits >
class FemWrapperTraits
{
public:
  typedef FemWrapper< SpaceTraits >                                   derived_type;
  typedef typename SpaceTraits::derived_type                          SpaceType;
  typedef typename SpaceTraits::BackendType::MapperType               BackendType;
  typedef typename SpaceTraits::GridPartType::IndexSetType::IndexType IndexType;
};


template< class SpaceTraits >
class FemWrapper
  : public Interface< FemWrapperTraits< SpaceTraits > >
{
public:
  typedef FemWrapper< SpaceTraits >       ThisType;
  typedef FemWrapperTraits< SpaceTraits > Traits;
  typedef Interface< Traits >             InterfaceType;

  typedef typename Traits::SpaceType    SpaceType;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::IndexType    IndexType;

  FemWrapper(const SpaceType& _space)
    : space_(_space)
  {}

  const SpaceType& space() const
  {
    return space_;
  }

  const BackendType& backend() const
  {
    return space_.backend().mapper();
  }

  template< class EntityType >
  IndexType numDofs(const EntityType& entity) const
  {
    return backend().numDofs(entity);
  }

  IndexType maxNumDofs() const
  {
    return backend().maxNumDofs();
  }

  class Functor
  {
  public:
    Functor(Dune::DynamicVector< IndexType >& globalIndices)
      : globalIndices_(globalIndices)
    {}

    void operator()(int localDoF,int globalDoF)
    {
      assert(localDoF < int(globalIndices_.size()));
      globalIndices_[localDoF] = globalDoF;
    }
  private:
    Dune::DynamicVector< IndexType >& globalIndices_;
  };

  template< class EntityType >
  void mapToGlobal(const EntityType& entity, Dune::DynamicVector< IndexType >& globalIndices) const
  {
    // some checks
    const IndexType numLocalDofs = numDofs(entity);
    if (globalIndices.size() < numLocalDofs)
      globalIndices.resize(numLocalDofs);
    // compute
    Functor functor(globalIndices);
    backend().mapEachEntityDof(entity, functor);
  }

  /**
   *  \attention  This method is implemented using mapToGlobal(entity, globalIndices) and thus not optimal!
   */
  template< class EntityType >
  IndexType mapToGlobal(const EntityType& entity, const IndexType& localIndex) const
  {
    const IndexType numLocalDofs = numDofs(entity);
    assert(localIndex < numLocalDofs);
    Dune::DynamicVector< IndexType > globalIndices(numLocalDofs);
    mapToGlobal(entity, globalIndices);
    return globalIndices[localIndex];
  }

private:
  const SpaceType& space_;
}; // class FemWrapper


} // namespace Mapper
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH
