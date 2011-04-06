#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH

namespace Dune
{

namespace Functionals
{

namespace LocalOperation
{

template< class FunctionSpaceImp, class Implementation >
class Interface
{
public:

  typedef Implementation
    ImplementationType;

  typedef FunctionSpaceImp
    FunctionSpaceType;

  enum
  {
    dimDomain = FunctionSpaceType::dimDomain,
    dimRange = FunctionSpaceType::dimRange
  };

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

  Interface()
  {
    std::cout << "Interface::Interface()" << std::endl;
  }

  ~Interface()
  {
    std::cout << "Interface::~Interface()" << std::endl;
  }

  template< class LocalFunctionType >
  RangeFieldType operateLocal( const LocalFunctionType& localFunction ) const
  {
    std::cout << "Interface::operateLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION( asImp().operateLocal( localFunction ) );
    return asImp().operateLocal( localFunction );
  }

  template< class LocalFunctionType, class LocalPointType >
  RangeFieldType evaluateLocal( const LocalFunctionType& localFunction,
                                const LocalPointType& localPoint ) const
  {
    std::cout << "Interface::evaluateLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION( asImp().evaluateLocal( localFunction, localPoint ) );
    return asImp().evaluateLocal( localFunction, localPoint );
  }

private:
    //! for CRTP trick
    ImplementationType& asImp()
    {
        return static_cast< ImplementationType& >(*this);
    }

    //! for CRTP trick
    const ImplementationType& asImp() const
    {
        return static_cast< const ImplementationType& >(*this);
    }

}; // end class Interface

} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
