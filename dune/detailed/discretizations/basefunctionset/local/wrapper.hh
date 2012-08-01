#ifndef DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_WRAPPER_HH
#define DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_WRAPPER_HH

namespace Dune
{

namespace Detailed {

namespace Discretizations
{

namespace BaseFunctionSet
{

namespace Local
{

template< class LocalFunctionImp >
class Wrapper
{
public:
  typedef LocalFunctionImp
    LocalFunctionType;

  typedef Wrapper< LocalFunctionType >
    ThisType;

  typedef typename LocalFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType
    BaseFunctionSetType;

  typedef typename DiscreteFunctionSpaceType::FunctionSpaceType
    FunctionSpaceType;

  typedef typename DiscreteFunctionSpaceType::GridPartType
    GridPartType;

  typedef typename DiscreteFunctionSpaceType::EntityType
    EntityType;

  enum{ polynomialOrder = DiscreteFunctionSpaceType::polynomialOrder };

  typedef typename FunctionSpaceType::DomainFieldType
    DomainFieldType;

  typedef typename FunctionSpaceType::DomainType
    DomainType;

  typedef typename FunctionSpaceType::RangeFieldType
    RangeFieldType;

  typedef typename FunctionSpaceType::RangeType
    RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType
    JacobianRangeType;

  typedef typename FunctionSpaceType::HessianRangeType
    HessianRangeType;

  //! constructor
  Wrapper( const LocalFunctionType& localFunction )
    : localFunction_( localFunction )
  {
  }

private:
  //! copy constructor
  Wrapper( const ThisType& );

public:
  const BaseFunctionSetType& baseFunctionSet() const
  {
    return localFunction_.space().baseFunctionSet();
  }

  const EntityType& entity() const
  {
    return localFunction_.entity();
  }

  unsigned int size() const
  {
    return 1;
  }

  int order() const
  {
    return localFunction_.order();
  }

  void evaluate( const DomainType& x, std::vector< RangeType >& ret) const
  {
    localFunction_.evaluate( x, ret[0] );
  }

  void jacobian( const DomainType& x, std::vector< JacobianRangeType >& ret ) const
  {
    localFunction_.jacobian( x, ret[0] );
  }

private:

  //! assignment operator
  ThisType& operator=( const ThisType& );

  const LocalFunctionType& localFunction_;

}; // end class Lagrange

} // end namespace Local

} // end namespace BaseFunctionSet

} // namespace Discretizations

} // namespace Detailed

} // end namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_BASEFUNCTIONSET_LOCAL_WRAPPER_HH
