#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

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

  ContinuousFiniteElement(GridPartType& gridPart)
    : gridPart_(gridPart)
    , hostSpace_(gridPart)
  {
  }

  const GridPartType& gridPart() const
  {
    return gridPart_;
  }

  const HostSpaceType& hostSpace() const
  {
    return hostSpace_;
  }

private:
  const GridPartType& gridPart_;
  const HostSpaceType hostSpace_;

}; // end class ContinuousFiniteElement

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
