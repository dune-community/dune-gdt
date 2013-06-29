#ifndef DUNE_GDT_BASEFUNCTIONSET_FEM_HH
#define DUNE_GDT_BASEFUNCTIONSET_FEM_HH

#ifdef HAVE_CMAKE_CONFIG
  #include "cmake_config.h"
#elif defined (HAVE_CONFIG_H)
  #include "config.h"
#endif

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace BaseFunctionSet {


// forward, to be used in the traits and to allow for specialization
template< class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemWrapper;


template< class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim, int rangeDimCols = 1 >
class FemWrapperTraits
{
public:
  typedef FemWrapper< FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, rangeDimCols > derived_type;
  typedef typename Dune::Fem::BaseFunctionSetInterface< FemBaseFunctionSetTraits >::BaseFunctionSetType BackendType;
  typedef EntityImp                                                                                     EntityType;
};


template< class FemBaseFunctionSetTraits, class EntityImp, class DomainFieldImp, int domainDim, class RangeFieldImp, int rangeDim >
class FemWrapper< FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
  : public BaseFunctionSetInterface< FemWrapperTraits< FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 >
{
public:
  typedef FemWrapperTraits< FemBaseFunctionSetTraits, EntityImp, DomainFieldImp, domainDim, RangeFieldImp, rangeDim, 1 > Traits;
  typedef typename Traits::BackendType  BackendType;
  typedef typename Traits::EntityType   EntityType;

  typedef DomainFieldImp                                  DomainFieldType;
  static const unsigned int                               dimDomain = domainDim;
  typedef Dune::FieldVector< DomainFieldType, dimDomain > DomainType;
  typedef RangeFieldImp                                 RangeFieldType;
  static const unsigned int                             dimRange = rangeDim;
  static const unsigned int                             dimRangeCols = 1;
  typedef Dune::FieldVector< RangeFieldType, dimRange > RangeType;
  typedef Dune::FieldMatrix< RangeFieldType, dimRange, dimDomain > JacobianRangeType;

  template< class S >
  FemWrapper(const Dune::Fem::DiscreteFunctionSpaceInterface< S >& femSpace, const EntityType& en)
    : entity_(en)
    , order_(femSpace.order())
    , backend_(femSpace.baseFunctionSet(entity_))
  {}

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

  void evaluate(const DomainType& x, std::vector< RangeType >& ret) const
  {
    assert(ret.size() >= backend_.size());
    backend_.evaluateAll(x, ret);
  }

  void jacobian(const DomainType& x, std::vector< JacobianRangeType >& ret) const
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
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_BASEFUNCTIONSET_FEM_HH
