#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LocalOperator {


template< class Traits >
class Codim0Interface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BinaryEvaluationType     BinaryEvaluationType;
  typedef typename Traits::LocalizableFunctionType  LocalizableFunctionType;

  const LocalizableFunctionType& inducingFunction() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().inducingFunction());
    return asImp().inducingFunction();
  }

  const BinaryEvaluationType& inducingEvaluation() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().inducingEvaluation());
    return asImp().inducingEvaluation();
  }

  size_t numTmpObjectsRequired() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numTmpObjectsRequired());
    return asImp().numTmpObjectsRequired();
  }

  /**
   *  \brief      Applies the local operator.
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A       Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the of the {testBase,ansatzBase}
   *  \tparam rC{T,a} dimRangeCols of the {testBase,ansatzBase}
   *  \attention  ret is assumed to be zero!
   */
  template< class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void apply(const BaseFunctionSetInterface< T, D, d, R, rT, rCT >& testBase,
             const BaseFunctionSetInterface< A, D, d, R, rA, rCA >& ansatzBase,
             Dune::DynamicMatrix< R >& ret,
             std::vector< Dune::DynamicMatrix< R > >& tmpLocalMatrices) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().apply(testBase, ansatzBase, ret, tmpLocalMatrices));
  }

  derived_type& asImp()
  {
    return static_cast< derived_type& >(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class Codim0Interface


} // namespace LocalOperator
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH
