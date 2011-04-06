#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH

namespace Dune {

namespace Functionals {

namespace LocalOperation {

template <class FunctionSpaceImp, class Implementation>
class Interface
{
public:
  typedef Implementation ImplementationType;

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

  Interface()
  {
    std::cout << "Interface::Interface()" << std::endl;
  }

  ~Interface()
  {
    std::cout << "Interface::~Interface()" << std::endl;
  }

  template <class LocalFunctionType>
  RangeFieldType operateLocal(const LocalFunctionType& localFunction) const
  {
    std::cout << "Interface::operateLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION(asImp().operateLocal(localFunction));
    return asImp().operateLocal(localFunction);
  }

  template <class LocalFunctionType, class LocalPointType>
  RangeFieldType evaluateLocal(const LocalFunctionType& localFunction, const LocalPointType& localPoint) const
  {
    std::cout << "Interface::evaluateLocal()" << std::endl;
    CHECK_INTERFACE_IMPLEMENTATION(asImp().evaluateLocal(localFunction, localPoint));
    return asImp().evaluateLocal(localFunction, localPoint);
  }

private:
  //! for CRTP trick
  ImplementationType& asImp()
  {
    return static_cast<ImplementationType&>(*this);
  }

  //! for CRTP trick
  const ImplementationType& asImp() const
  {
    return static_cast<const ImplementationType&>(*this);
  }

}; // end class Interface


template <class Implementation>
class TypeSelector
{
public:
  typedef Interface<typename Implementation::FunctionSpaceType, Implementation> Select;

}; // end class InterfaceSelector


template <class Implementation>
class Wrapper : public Interface<typename Implementation::FunctionSpaceType, Wrapper<Implementation>>
{
public:
  typedef Interface<typename Implementation::FunctionSpaceType, Wrapper<Implementation>> InterfaceType;

  typedef Implementation ImplementationType;

  typedef typename ImplementationType::FunctionSpaceType FunctionSpaceType;

  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  Wrapper(const ImplementationType& implementation)
    : InterfaceType()
    , implementation_(implementation)
  {
  }

  template <class LocalFunctionType>
  RangeFieldType operateLocal(const LocalFunctionType& localFunction) const
  {
    return implementation_.operateLocal(localFunction);
  }

  template <class LocalFunctionType, class LocalPointType>
  RangeFieldType evaluateLocal(const LocalFunctionType& localFunction, const LocalPointType& localPoint) const
  {
    return implementation_.evaluateLocal(localFunction, localPoint);
  }

private:
  const ImplementationType implementation_;
}; // end class

} // end namespace LocalOperation

} // end namespace Functionals

} // end namespace Dune

#endif // DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
