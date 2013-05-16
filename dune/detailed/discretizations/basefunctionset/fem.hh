#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_HH

#ifdef HAVE_CMAKE_CONFIG
#include "cmake_config.h"
#elif defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template <class FemBaseFunctionSetTraits, class EntityImp>
class FemWrapper;


template <class FemBaseFunctionSetTraits, class EntityImp>
class FemWrapperTraits
{
public:
  typedef FemWrapper<FemBaseFunctionSetTraits, EntityImp> derived_type;
  typedef typename Dune::Fem::BaseFunctionSetInterface<FemBaseFunctionSetTraits>::BaseFunctionSetType BackendType;
  typedef EntityImp EntityType;
  typedef typename BackendType::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = BackendType::dimDomain;
  typedef typename BackendType::DomainType DomainType;
  typedef typename BackendType::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = BackendType::dimRange;
  typedef typename BackendType::RangeType RangeType;
  typedef typename BackendType::JacobianRangeType JacobianRangeType;
};


template <class FemBaseFunctionSetTraits, class EntityImp>
class FemWrapper : public BaseFunctionSetInterface<FemWrapperTraits<FemBaseFunctionSetTraits, EntityImp>>
{
public:
  typedef FemWrapperTraits<FemBaseFunctionSetTraits, EntityImp> Traits;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  typedef typename Traits::DomainFieldType DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  static const unsigned int dimRange = Traits::dimRange;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::JacobianRangeType JacobianRangeType;

  template <class S>
  FemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface<S>& femSpace, const EntityType& en)
    : entity_(en)
    , order_(femSpace.order())
    , backend_(femSpace.baseFunctionSet(entity_))
  {
  }

  const EntityType& entity() const
  {
    return entity_;
  }

  const BackendType& backend() const
  {
    return backend_;
  }

  size_t size() const
  {
    return backend_.size();
  }

  size_t order() const
  {
    return order_;
  }

  void evaluate(const DomainType& x, std::vector<RangeType>& ret) const
  {
    assert(ret.size() >= backend_.size());
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector<JacobianRangeType>& ret) const
  {
    assert(ret.size() >= backend_.size());
    backend_.jacobianAll(x, entity_.geometry().jacobianInverseTransposed(x), ret);
  }

private:
  const EntityType& entity_;
  const size_t order_;
  const BackendType backend_;
}; // class FemWrapper


} // namespace BaseFunctionSet
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_FEM_HH
