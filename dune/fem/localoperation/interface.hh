#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH

namespace Dune
{

namespace Functionals
{

namespace LocalOperation
{

template< class FunctionSpaceImp >
class Interface
{
public:

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
  RangeFieldType operate( const LocalFunctionType& localFunction ) const
  {
    std::cout << "Interface::operateLocal()" << std::endl;
  }

  template< class LocalFunctionType, class LocalPointType >
  RangeFieldType evaluate( const LocalFunctionType& localFunction,
                                const LocalPointType& localPoint ) const
  {
    std::cout << "Interface::evaluateLocal()" << std::endl;
  }
}; // end class Interface


} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
