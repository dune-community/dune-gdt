#ifndef DUNE_DETAILED_DISCRETIZATIONS_SPACE_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_SPACE_INTERFACE_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


template< class Traits >
class SpaceInterface
{
public:
  typedef SpaceInterface< Traits >      ThisType;
  typedef typename Traits::derived_type derived_type;

  typedef typename Traits::GridPartType   GridPartType;
  typedef typename GridPartType::ctype    DomainFieldType;
  static const unsigned int               dimDomain = GridPartType::dimension;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int               dimRange = Traits::dimRange;
  static const int                        polOrder = Traits::polynomialOrder;

  typedef typename Traits::BackendType          BackendType;
  typedef typename Traits::MapperType           MapperType;
  typedef typename Traits::BaseFunctionSetType  BaseFunctionSetType;
  typedef typename Traits::EntityType           EntityType;

  const GridPartType& gridPart() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().gridPart());
    return asImp().gridPart();
  }

  const BackendType& backend() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().backend());
    return asImp().backend();
  }

  bool continuous() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().continuous());
    return asImp().continuous();
  }

  const MapperType& mapper() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().mapper());
    return asImp().mapper();
  }

  BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().baseFunctionSet(entity));
    return asImp().baseFunctionSet(entity);
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class SpaceInterface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_SPACE_INTERFACE_HH
