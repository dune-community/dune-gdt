#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


template <class Traits>
class LocalOperatorCodim0Interface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BinaryEvaluationType BinaryEvaluationType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  const LocalizableFunctionType& inducingFunction() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().inducingFunction());
    return asImp().inducingFunction();
  }

  size_t numTmpObjectsRequired() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numTmpObjectsRequired());
    return asImp().numTmpObjectsRequired();
  }

  /**
   *  \brief      Applies the local operator.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \tparam A   Traits of the BaseFunctionSetInterface implementation, representing the type of the ansatzBase
   *  \attention  ret is assumed to be zero!
   */
  template <class T, class A, class RangeFieldType>
  void apply(const BaseFunctionSetInterface<T>& testBase, const BaseFunctionSetInterface<A>& ansatzBase,
             Dune::DynamicMatrix<RangeFieldType>& ret,
             std::vector<Dune::DynamicMatrix<RangeFieldType>>& tmpLocalMatrices) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().apply(testBase, ansatzBase, ret, tmpLocalMatrices));
    asImp().apply(testBase, ansatzBase, ret, tmpLocalMatrices);
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class LocalOperatorCodim0Interface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALOPERATOR_INTERFACE_HH
