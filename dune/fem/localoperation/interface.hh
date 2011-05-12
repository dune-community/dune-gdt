#ifndef DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH
#define DUNE_FEM_FUNCTIONALS_LOCALOPERATION_INTERFACE_HH

namespace Dune {

namespace Functionals {

namespace LocalOperation {

/**
 * @brief Interface class for local operations.
 *
 * Let @f$V_1@f$, @f$V_2@f$ be local function spaces and
 * @f$v_1\in V_1@f$, @f$v_2\in V_2@f$ local functions.
 * A @e local @e operation @f$L@f$ can be either defined by
 * @f$L:V_1\rightarrow \mathbb{R}@f$
 * or @f$L:V_1\times V_2\rightarrow\mathbb{R}@f$.
 *
 * @tparam FunctionSpaceImp The function space.
 */
template <class FunctionSpaceImp>
class Interface
{
public:
  //! Type of the function space.
  typedef FunctionSpaceImp FunctionSpaceType;

  enum
  {
    //! Dimension of the domain.
    dimDomain = FunctionSpaceType::dimDomain,
    //! Dimension of the range.
    dimRange = FunctionSpaceType::dimRange
  };

  //! Intrinsic type used for values in the domain field (usually a double).
  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;

  //! Type of domain vector (using type of domain field) has a Dune::FieldVector type interface.
  typedef typename FunctionSpaceType::DomainType DomainType;

  //! Intrinsic type used for values in the range field (usually a double).
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  //! Type of domain vector (using type of range field) has a Dune::FieldVector type interface.
  typedef typename FunctionSpaceType::RangeType RangeType;

  //! Intrinsic type used for the jacobian values has a Dune::FieldMatrix type interface.
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

  /**
   * @brief Interface method for ???.
   *
   * @param localTestFunction
   *
   * @todo Please doc me!
   */
  template <class LocalTestFunctionType>
  RangeFieldType operate(const LocalTestFunctionType& localTestFunction) const
  {
    std::cout << "LocalOperation::Interface::operate()" << std::endl;
  }

  /**
   * @brief Interface method for evaluating the local operation.
   *
   * @param localTestFunction The local test function.
   * @param localPoint The local point, where we want to evaluate the local operation.
   */
  template <class LocalTestFunctionType, class LocalPointType>
  RangeFieldType evaluate(const LocalTestFunctionType& localTestFunction, const LocalPointType& localPoint) const
  {
    std::cout << "LocalOperation::Interface::evaluate()" << std::endl;
  }

  /**
   * @brief Interface method for ???.
   *
   * @param localAnsatzFunction The local ansatz function.
   * @param localTestFunction The local test function.
   *
   * @todo Please doc me!
   */
  template <class LocalAnsatzFunctionType, class LocalTestFunctionType>
  RangeFieldType operate(const LocalAnsatzFunctionType& localAnsatzFunction,
                         const LocalTestFunctionType& localTestFunction) const
  {
    std::cout << "LocalOperation::Interface::operate()" << std::endl;
  }

  /**
   * @brief Interface method for evaluating the local operation.
   *
   * @param localAnsatzFunction The local ansatz function.
   * @param localTestFunction The local test function.
   * @param localPoint The local point, where we want to evaluate the local operation.
   */
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
