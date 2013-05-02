#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH

#include <dune/common/dynvector.hh>

#include <dune/fem/space/mapper/dofmapper.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {


//// forward, to be used in the traits and to allow for specialization
template< class FemDofMapperImp >
class MapperWrappedFemDofMapper;


template< class FemDofMapperImp >
class MapperWrappedFemDofMapperTraits
{
public:
  typedef MapperWrappedFemDofMapper< FemDofMapperImp >              derived_type;
  typedef Dune::Fem::DofMapper< typename FemDofMapperImp::Traits >  BackendType;
};


template< class FemDofMapperImp >
class MapperWrappedFemDofMapper
  : public MapperInterface< MapperWrappedFemDofMapperTraits< FemDofMapperImp > >
{
public:
  typedef MapperWrappedFemDofMapperTraits< FemDofMapperImp > Traits;
  typedef typename Traits::BackendType    BackendType;

  MapperWrappedFemDofMapper(const BackendType& femMapper)
    : backend_(femMapper)
  {}

  const BackendType& backend() const
  {
    return backend_;
  }

  template< class EntityType >
  size_t numDofs(const EntityType& entity) const
  {
    return backend_.numDofs(entity);
  }

  size_t maxNumDofs() const
  {
    return backend_.maxNumDofs();
  }

private:
  class Functor
  {
  public:
    Functor(Dune::DynamicVector< size_t >& globalIndices)
      : globalIndices_(globalIndices)
    {}

    void operator()(int localDoF, int globalDoF)
    {
      assert(localDoF < globalIndices_.size());
      globalIndices_[localDoF] = globalDoF;
    }
  private:
    Dune::DynamicVector< size_t >& globalIndices_;
  };

public:
  template< class EntityType >
  void mapToGlobal(const EntityType& entity, Dune::DynamicVector< size_t >& globalIndices) const
  {
    // some checks
    const size_t numLocalDofs = numDofs(entity);
    if (globalIndices.size() < numLocalDofs)
      globalIndices.resize(numLocalDofs);
    // compute
    Functor functor(globalIndices);
    backend_.mapEachEntityDof(entity, functor);
  }

  /**
   *  \attention  This method is implemented using mapToGlobal(entity, globalIndices) and thus not optimal!
   */
  template< class EntityType >
  size_t mapToGlobal(const EntityType& entity, const size_t& localIndex) const
  {
    const size_t numLocalDofs = numDofs(entity);
    assert(localIndex < numLocalDofs);
    Dune::DynamicVector< size_t > globalIndices(numLocalDofs);
    mapToGlobal(entity, globalIndices);
    return globalIndices[localIndex];
  }

private:
  const BackendType& backend_;
}; // class MapperWrappedFemDofMapper


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_FEM_LOCALFUNCTIONS_HH
