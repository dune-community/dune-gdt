#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH

#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


template< class Traits >
class BinaryEvaluationInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes a binary evaluation.
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A Traits of the ansatz BaseFunctionSetInterface implementation
   */
  template< class L, class T, class A, class D, int d, class R, int rR, int rC >
  static void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rR, rC >& localFunction,
                       const BaseFunctionSetInterface< T >& testBase,
                       const BaseFunctionSetInterface< A >& ansatzBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicMatrix< R >& ret)
  {
    derived_type::evaluate(localFunction, testBase, ansatzBase, localPoint, ret);
  }
}; // class BinaryEvaluationInterface


template< class Traits >
class UnaryEvaluationInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes a unary evaluation.
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   */
  template< class L, class T, class D, int d, class R, int rR, int rC >
  static void evaluate(const Dune::Stuff::LocalFunctionInterface< L, D, d, R, rR, rC >& localFunction,
                       const BaseFunctionSetInterface< T >& testBase,
                       const Dune::FieldVector< D, d >& localPoint,
                       Dune::DynamicVector< R >& ret)
  {
    derived_type::evaluate(localFunction, testBase, localPoint, ret);
  }
}; // class UnaryEvaluationInterface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH
