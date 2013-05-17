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


template <class Traits>
class UnaryEvaluationInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes a unary evaluation.
   *  \tparam L       Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{L,T}  dimRange of the {localfunction,testBase}
   *  \tparam rC{L,T} dimRangeRows of the {localfunction,testBase}
   */
  template <class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
                       const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret)
  {
    derived_type::evaluate(localFunction, testBase, localPoint, ret);
  }
}; // class UnaryEvaluationInterface


template <class Traits>
class BinaryEvaluationInterface
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes a binary evaluation.
   *  \tparam L         Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T         Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A         Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localfunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localfunction,testBase,ansatzBase}
   */
  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
                       const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                       const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<R>& ret)
  {
    derived_type::evaluate(localFunction, testBase, ansatzBase, localPoint, ret);
  }
}; // class BinaryEvaluationInterface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH
