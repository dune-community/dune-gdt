#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH

namespace Dune {

namespace Functionals {

namespace LocalOperation {

template <class FunctionSpaceImp>
class Interface
{
public:
  typedef FunctionSpaceImp FunctionSpaceType;

  enum
  {
    dimDomain = FunctionSpaceType::dimDomain,
    dimRange  = FunctionSpaceType::dimRange
  };

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::RangeType RangeType;

  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  template <class LocalTestFunctionType>
  RangeFieldType operate(const LocalTestFunctionType& localTestFunction) const
  {
    std::cout << "LocalOperation::Interface::operate()" << std::endl;
  }

  template <class LocalTestFunctionType, class LocalPointType>
  RangeFieldType evaluate(const LocalTestFunctionType& localTestFunction, const LocalPointType& localPoint) const
  {
    std::cout << "LocalOperation::Interface::evaluate()" << std::endl;
  }

  template <class LocalAnsatzFunctionType, class LocalTestFunctionType>
  RangeFieldType operate(const LocalAnsatzFunctionType& localAnsatzFunction,
                         const LocalTestFunctionType& localTestFunction) const
  {
    std::cout << "LocalOperation::Interface::operate()" << std::endl;
  }

  template <class LocalAnsatzFunctionType, class LocalTestFunctionType, class LocalPointType>
  RangeFieldType evaluate(const LocalAnsatzFunctionType& localAnsatzFunction,
                          const LocalTestFunctionType& localTestFunction, const LocalPointType& localPoint) const
  {
    std::cout << "LocalOperation::Interface::evaluate()" << std::endl;
  }

}; // end class Interface


} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
