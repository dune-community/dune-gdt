#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

// dune-functionals includes
#include <dune/functionals/common/localbasefunctionprovider.hh>

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

  typedef typename HostSpaceType::RangeFieldType RangeFieldType;

  typedef typename HostSpaceType::EntityType EntityType;

  typedef Dune::Functionals::Common::LocalBaseFunctionProvider<HostSpaceType> LocalBaseFunctionProviderType;

  ContinuousFiniteElement(GridPartType& gridPart)
    : gridPart_(gridPart)
    , hostSpace_(gridPart)
    , localBaseFunctionProvider_(hostSpace_)
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

  const LocalBaseFunctionProviderType& localBaseFunctionProvider() const
  {
    return localBaseFunctionProvider_;
  }

private:
  const GridPartType& gridPart_;
  const HostSpaceType hostSpace_;
  const LocalBaseFunctionProviderType localBaseFunctionProvider_;

}; // end class ContinuousFiniteElement

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_FINITEELEMENT_HH
