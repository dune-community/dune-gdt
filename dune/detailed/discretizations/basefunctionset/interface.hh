#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace BaseFunctionSet {


template <class Traits>
class Interface
{
public:
  typedef Interface<Traits> ThisType;
  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::IndexType IndexType;
  typedef typename Traits::EntityType EntityType;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::JacobianRangeType JacobianRangeType;

  const SpaceType& space() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().space());
    return asImp().space();
  }

  const EntityType& entity() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().entity());
    return asImp().entity();
  }

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  IndexType size() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().size());
    return asImp().size();
  }

  unsigned int order() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order());
    return asImp().order();
  }

  void evaluate(const DomainType& x, std::vector<RangeType>& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().evaluate(x, ret));
    return asImp().evaluate(x, ret);
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().jacobian(x, ret));
    return asImp().jacobian(x, ret);
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


} // namespace BaseFunctionSet
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_INTERFACE_HH
