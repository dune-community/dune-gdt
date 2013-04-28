#ifndef DUNE_DETAILED_DISCRETIZATIONS_MAPPER_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_MAPPER_INTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynvector.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace Mapper {


template <class Traits>
class Interface
{
public:
  typedef Interface<Traits> ThisType;
  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::IndexType IndexType;

  const SpaceType& space() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().space());
    return asImp().space();
  }

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  IndexType maxNumDofs() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().maxNumDofs());
    return asImp().maxNumDofs();
  }

  template <class EntityType>
  IndexType numDofs(const EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numDofs(entity));
    return asImp().numDofs(entity);
  }

  template <class EntityType>
  void mapToGlobal(const EntityType& entity, Dune::DynamicVector<IndexType>& globalIndices) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity, globalIndices));
    asImp().mapToGlobal(entity, globalIndices);
  }

  template <class EntityType>
  IndexType mapToGlobal(const EntityType& entity, const IndexType& localIndex) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapToGlobal(entity, localIndex));
    return asImp().mapToGlobal(entity, localIndex);
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Interface


} // namespace Mapper
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune


#endif // DUNE_DETAILED_DISCRETIZATIONS_MAPPER_INTERFACE_HH
