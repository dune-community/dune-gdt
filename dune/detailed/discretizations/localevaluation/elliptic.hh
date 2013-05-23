#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_ELLIPTIC_HH

#include <tuple>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/static_assert.hh>

#include <dune/stuff/function/interface.hh>
#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

#include "interface.hh"

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LocalEvaluation {


// forward, to be used in the traits
template <class LocalizableFunctionImp>
class Elliptic;


/**
 *  \brief  Traits for the Elliptic evaluation.
 */
template <class LocalizableFunctionImp>
class EllipticTraits
{
public:
  typedef Elliptic<LocalizableFunctionImp> derived_type;
  typedef LocalizableFunctionImp LocalizableFunctionType;
  dune_static_assert((Dune::IsBaseOf<Dune::Stuff::LocalizableFunction, LocalizableFunctionImp>::value),
                     "ERROR: LocalizableFunctionImp is not a Dune::Stuff::LocalizableFunction.");
};


/**
 *  \brief  Computes an elliptic evaluation.
 */
template <class LocalizableFunctionImp>
class Elliptic : public LocalEvaluation::Codim0Interface<EllipticTraits<LocalizableFunctionImp>, 2>
{
public:
  typedef EllipticTraits<LocalizableFunctionImp> Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;

  Elliptic(const LocalizableFunctionType& inducingFunction)
    : inducingFunction_(inducingFunction)
  {
  }

  //  template< class EntityType >
  //  auto localFunctions(const EntityType& entity) const -> decltype(LocalizableFunctionType::template localFunction<
  //  EntityType >(entity))
  //  {
  //    return inducingFunction_.localFunction(entity);
  //  }

  template <class EntityType>
  std::tuple<typename LocalizableFunctionType::template LocalFunction<EntityType>::Type>
  localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(inducingFunction_.localFunction(entity));
  }

  /**
   * \brief extracts the local functions and calls the correct order() method
   */
  template <class L, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  int order(const std::tuple<L>& localFunctions, const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    const auto& localFunction = std::get<0>(localFunctions);
    return order(localFunction, testBase, ansatzBase);
  }

  /**
   *  \todo add copydoc
   *  \return localFunction.order() + (testBase.order() - 1) + (ansatzBase.order() - 1)
   */
  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  int order(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
            const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    if (localFunction.order() < 0)
      return -1;
    else
      return std::max(int(localFunction.order() + testBase.order() + ansatzBase.order() - 2), 0);
  } // int order(...)

  /**
   * \brief extracts the local functions and calls the correct evaluate() method
   */
  template <class L, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  static void evaluate(const std::tuple<L>& localFunctions,
                       const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                       const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<R>& ret)
  {
    const auto& localFunction = std::get<0>(localFunctions);
    evaluate(localFunction, testBase, ansatzBase, localPoint, ret);
  }

  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& /*localFunction*/,
                       const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& /*testBase*/,
                       const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& /*ansatzBase*/,
                       const Dune::FieldVector<D, d>& /*localPoint*/, Dune::DynamicMatrix<R>& /*ret*/)
  {
    dune_static_assert((Dune::AlwaysFalse<R>::value), "ERROR: not implemented for this combination of dimensions!");
  }

  /**
   *  \brief  Computes an elliptic evaluationf for a scalar local function and scalar basefunctionsets.
   *  \tparam L Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D DomainFieldType
   *  \tparam d dimDomain
   *  \tparam R RangeFieldType
   */
  template <class L, class T, class A, class D, int d, class R>
  static void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, 1, 1>& localFunction,
                       const BaseFunctionSetInterface<T, D, d, R, 1, 1>& testBase,
                       const BaseFunctionSetInterface<A, D, d, R, 1, 1>& ansatzBase,
                       const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<R>& ret)
  {
    typedef typename BaseFunctionSetInterface<A, D, d, R, 1, 1>::RangeType RangeType;
    typedef typename BaseFunctionSetInterface<A, D, d, R, 1, 1>::JacobianRangeType JacobianRangeType;
    // evaluate local function
    const RangeType functionValue = localFunction.evaluate(localPoint);
    // evaluate test gradient
    const size_t rows = testBase.size();
    std::vector<JacobianRangeType> testGradients(rows, JacobianRangeType(0));
    testBase.jacobian(localPoint, testGradients);
    // evaluate ansatz gradient
    const size_t cols = ansatzBase.size();
    std::vector<JacobianRangeType> ansatzGradients(cols, JacobianRangeType(0));
    ansatzBase.jacobian(localPoint, ansatzGradients);
    // compute products
    assert(ret.rows() >= rows);
    assert(ret.cols() >= cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      auto& retRow = ret[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const R gradientProduct = ansatzGradients[jj][0] * testGradients[ii][0];
        retRow[jj]              = functionValue * gradientProduct;
      }
    }
  } // ... evaluate< ..., 1, 1 >(...)

private:
  const LocalizableFunctionType& inducingFunction_;
}; // class LocalElliptic


} // namespace Evaluation
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif //
