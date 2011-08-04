#ifndef DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGEFEMWRAPPER_HH
#define DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGEFEMWRAPPER_HH

// dune-fem includes
#include <dune/fem/space/lagrangespace.hh>

// dune-functionals include
#include <dune/functionals/discretefunctionspace/continuous/lagrange.hh>

namespace Dune
{

namespace Functionals
{

namespace DiscreteFunctionSpace
{

namespace Continuous
{

template< class HostSpaceImp >
class LagrangeFemAdapter
  : public DiscreteFunctionSpaceDefault<  LagrangeDiscreteFunctionSpaceTraits<  typename HostSpaceImp::FunctionSpaceType,
                                                                                typename HostSpaceImp::GridPartType,
                                                                                HostSpaceImp::polynomialOrder,
                                                                                CachingStorage > >
{
public:
  //! this should be of type Dune::Functionals::DiscreteFunctionSpace::Continuous::Lagrange
  typedef HostSpaceImp
    HostSpaceType;

  typedef typename HostSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename HostSpaceType::GridPartType
    GridPartType;

  enum{ polOrder = HostSpaceType::polynomialOrder };

  enum{ polynomialOrder = HostSpaceType::polynomialOrder };

  typedef LagrangeDiscreteFunctionSpaceTraits<  FunctionSpaceType,
                                                GridPartType,
                                                polOrder,
                                                CachingStorage >
    Traits;

  typedef LagrangeFemAdapter< HostSpaceType >
    LagrangeDiscreteFunctionSpaceType;

  typedef typename GridPartType::GridType
    GridType;

  typedef typename Traits::IndexSetType
    IndexSetType;

  typedef typename Traits::IteratorType
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

  typedef typename HostSpaceType::LocalBaseFunctionSetType
    HostLocalBaseFunctionSetType;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionSpaceType
    BaseFunctionSpaceType;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionSetImp
    BaseFunctionSetImp;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionMapType
    BaseFunctionMapType;

  typedef typename HostLocalBaseFunctionSetType::ScalarFactoryType
    ScalarFactoryType;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionSetSingletonFactoryType
    BaseFunctionSetSingletonFactoryType;

  typedef typename HostLocalBaseFunctionSetType::BaseFunctionSetSingletonProviderType
    BaseFunctionSetSingletonProviderType;

  typedef typename HostSpaceType::MapperType
    HostMapperType;

  typedef typename HostMapperType::LagrangePointSetType
    LagrangePointSetType;

  typedef typename HostMapperType::LagrangePointSetMapType
    LagrangePointSetMapType;

  typedef typename HostMapperType::MapperType
    MapperType;

  typedef typename HostMapperType::BlockMapperType
    BlockMapperType;

  enum{ localBlockSize = Traits::localBlockSize };

  typedef RangeFieldType
    DofType;

  enum{ dimVal = 1 };

  typedef DofManager< GridType >
    DofManagerType;

  typedef typename HostMapperType::MapperSingletonKeyType
    MapperSingletonKeyType;

  typedef LagrangeMapperSingletonFactory< MapperSingletonKeyType, MapperType >
    MapperSingletonFactoryType;

  typedef typename HostMapperType::BlockMapperSingletonFactoryType
    BlockMapperSingletonFactoryType;

  typedef typename HostMapperType::BlockMapperProviderType
    BlockMapperProviderType;

  typedef int
    IdentifierType;

private:
  typedef LagrangeDiscreteFunctionSpaceType
    ThisType;

  typedef Dune::DiscreteFunctionSpaceDefault< Traits >
    BaseType;

public:
  static const IdentifierType id = 665;

  using BaseType::gridPart;

  static const InterfaceType defaultInterface = InteriorBorder_InteriorBorder_Interface;

  static const CommunicationDirection defaultDirection = ForwardCommunication;

  explicit LagrangeFemAdapter(  const HostSpaceType& hostSpace,
                                const InterfaceType commInterface = defaultInterface,
                                const CommunicationDirection commDirection = defaultDirection )
    : BaseType( const_cast< GridPartType& >( hostSpace.gridPart() ), commInterface, commDirection ),
      hostSpace_( hostSpace )
  {
  }

  inline bool contains( const int codim ) const
  {
    return blockMapper().contains( codim );
  }

  inline bool continuous() const
  {
    return ( polynomialOrder > 0 );
  }

  inline DFSpaceIdentifier type() const
  {
    return LagrangeSpace_id;
  }

  inline int order() const
  {
    return polynomialOrder;
  }

  template< class EntityType >
  inline const BaseFunctionSetType baseFunctionSet( const EntityType& entity ) const
  {
    const HostLocalBaseFunctionSetType hostLocalBaseFunctionSet( hostSpace_, entity );
    return hostLocalBaseFunctionSet.baseFunctionSet_[entity.type()];
  }

  template< class EntityType >
  inline const LagrangePointSetType& lagrangePointSet( const EntityType& entity ) const
  {
    return hostSpace_.mapper_.lagrangePointSet_[entity.type()];
  }

  inline int dimensionOfValue() const
  {
    return dimVal;
  }

  MapperType& mapper() const
  {
    return *( hostSpace_.mapper_.mapper_ );
  }

  BlockMapperType& blockMapper() const
  {
    return *( hostSpace_.mapper_.blockMapper_ );
  }

private:

  //! copy constructor
  LagrangeFemAdapter( const ThisType& );

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const HostSpaceType& hostSpace_;

}; // end class LagrangeFemAdapter

} // end namespace Continuous

} // end namespace DiscreteFunctionSpace

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FUNCTIONALS_DISCRETEFUNCTIONSPACE_CONTINUOUS_LAGRANGEFEMWRAPPER_HH
