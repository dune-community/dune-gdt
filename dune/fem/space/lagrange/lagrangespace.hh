#ifndef DUNE_FEM_SPACE_LAGRANGESPACE_LAGRANGESPACE_HH
#define DUNE_FEM_SPACE_LAGRANGESPACE_LAGRANGESPACE_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace/lagrangespace.hh>

// dune-functionals includes
#include <dune/functionals/discretefunctionspace/continuous/lagrange.hh>

namespace Dune
{

namespace Fem
{

namespace Space
{

template< class FunctionalsLagrangeSpaceImp >
class LagrangeDiscreteFunctionSpace
//  : public DiscreteFunctionSpaceDefault<  LagrangeDiscreteFunctionSpaceTraits<  typename FunctionalsLagrangeSpaceImp::FunctionSpaceType,
//                                                                                typename FunctionalsLagrangeSpaceImp::GridPartType,
//                                                                                FunctionalsLagrangeSpaceImp::polynomialOrder,
//                                                                                CachingStorage > >
{
public:
  typedef FunctionalsLagrangeSpaceImp
    HostSpaceType;

  typedef typename HostSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename HostSpaceType::GridPartType
    GridPartType;

  enum{ polynomialOrder = HostSpaceType::polynomialOrder };

  typedef LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                GridPartType,
                                                polynomialOrder,
                                                CachingStorage >
    Traits;

  typedef LagrangeDiscreteFunctionSpace< HostSpaceType >
    LagrangeDiscreteFunctionSpaceType;

  typedef typename GridPartType::GridType
    GridType;

  typedef Traits::IndexSetType
    IndexSetType;

  typedef Traits::IteratorType
    IteratorType;

  enum{ dimension = GridType::dimension };

  typedef typename HostSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename HostSpaceType::DomainType
    DomainType;

  typedef typename HostSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename HostSpaceType::RangeType
    RangeType;

  enum{ dimRange = FunctionSpaceType::dimRange };

  typedef typename Traits::BaseFunctionSpaceType
    BaseFunctionSpaceType;

}; // end class LagrangeDiscreteFunctionSpace

} // end namespace Space

} // end namespace Fem

} // end namespace Dune


#endif // DUNE_FEM_SPACE_LAGRANGESPACE_LAGRANGESPACE_HH
