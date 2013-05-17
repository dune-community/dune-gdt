#ifndef DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynvector.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {


template <class Traits>
class LocalFunctionalCodim0Interface
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::UnaryEvaluationType UnaryEvaluationType;
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
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the BaseFunctionSetInterface implementation, representing the type of the testBase
   *  \attention  ret is assumed to be zero!
   */
  template <class T, class D, int d, class R, int r, int rC>
  void apply(const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase, Dune::DynamicVector<R>& ret,
             std::vector<Dune::DynamicVector<R>>& tmpLocalVectors) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().apply(testBase, ret, tmpLocalVectors));
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class LocalFunctionalCodim0Interface


} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_LOCALFUNCTIONAL_INTERFACE_HH
