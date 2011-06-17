#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

// dune-functionals includes
#include <dune/functionals/common/localbasefunction.hh>

namespace Dune {

namespace Functionals {

namespace DiscreteFunctionSpace {

template <class FunctionSpaceImp, class GridPartImp, int polOrder>
class ContinuousFiniteElement
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  typedef GridPartImp GridPartType;

  typedef Dune::LagrangeDiscreteFunctionSpace<FunctionSpaceType, GridPartType, polOrder> HostSpaceType;

  typedef Dune::Functionals::Common::LocalBaseFunction<HostSpaceType> LocalBaseFunctionType;

  typedef typename HostSpaceType::EntityType EntityType;

  typedef typename HostSpaceType::RangeFieldType RangeFieldType;

  typedef typename HostSpaceType::DomainType DomainType;

  /**
    \defgroup dune-fem related
    \{
    **/
  typedef typename HostSpaceType::BaseFunctionSetType BaseFunctionSetType;

  typedef typename HostSpaceType::IteratorType IteratorType;
  /**
    \}
    **/

  ContinuousFiniteElement(GridPartType& gridPart)
    : gridPart_(gridPart)
    , hostSpace_(gridPart)
    , numMaxDoFs_(-1)
  {
    // in the simple case, there should be the same number of dofs on each entity
    const IteratorType entityIterator         = begin();
    const EntityType& entity                  = *entityIterator;
    const BaseFunctionSetType baseFunctionSet = hostSpace_.baseFunctionSet(entity);
    numMaxDoFs_                               = baseFunctionSet.numBaseFunctions();
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const HostSpaceType& hostSpace() const
  {
    return hostSpace_;
  }

  const LocalBaseFunctionType localBaseFunction(const EntityType& entity, const unsigned int localDoFNumber) const
  {
    return LocalBaseFunctionType(entity, hostSpace_.baseFunctionSet(entity), localDoFNumber);
  }

  const unsigned int size() const
  {
    return hostSpace_.size();
  }

  const int numMaxDoFs() const
  {
    return numMaxDoFs_;
  }

  const int order() const
  {
    return hostSpace_.order();
  }

  /**
    \defgroup dune-fem related
    \{
    **/
  IteratorType begin() const
  {
    return hostSpace_.begin();
  }

  const IteratorType end() const
  {
    return hostSpace_.end();
  }

  const BaseFunctionSetType baseFunctionSet(const EntityType& entity) const
  {
    return hostSpace_.baseFunctionSet(entity);
  }

  const int mapToGlobal(const EntityType& entity, const int localDof) const
  {
    return hostSpace_.mapToGlobal(entity, localDof);
  }
  /**
    \}
    **/

private:
  const GridPartType& gridPart_;
  const HostSpaceType hostSpace_;
  unsigned int numMaxDoFs_;

}; // end class ContinuousFiniteElement

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
