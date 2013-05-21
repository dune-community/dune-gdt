#ifndef DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH
#define DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH

#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/stuff/localfunction/interface.hh>

#include <dune/detailed/discretizations/basefunctionset/interface.hh>

namespace Dune {
namespace Detailed {
namespace Discretizations {
namespace LocalEvaluation {


/**
 *  \brief  Interface for local evaluations that depend on a codim 0 entity.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 */
template <class Traits, int numArguments>
class Codim0Interface
{
public:
  Codim0Interface() = delete;
};


/**
 *  \brief  Interface for unary codim 0 evaluations.
 */
template <class Traits>
class Codim0Interface<Traits, 1>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam L       Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{L,T}  dimRange of the {localfunction,testBase}
   *  \tparam rC{L,T} dimRangeRows of the {localfunction,testBase}
   */
  template <class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  int order(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
            const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunction, testBase));
    asImp().order(localFunction, testBase);
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam L       Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{L,T}  dimRange of the {localfunction,testBase}
   *  \tparam rC{L,T} dimRangeRows of the {localfunction,testBase}
   *  \attention ret is assumed to be zero!
   */
  template <class L, class T, class D, int d, class R, int rL, int rCL, int rT, int rCT>
  void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
                const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunction, testBase, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Codim0Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 0 evaluations.
 */
template <class Traits>
class Codim0Interface<Traits, 2>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam L       Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{L,T}  dimRange of the {localfunction,testBase}
   *  \tparam rC{L,T} dimRangeRows of the {localfunction,testBase}
   */
  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  int order(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
            const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunction, testBase, ansatzBase));
    asImp().order(localFunction, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 0 evaluation.
   *  \tparam L         Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T         Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A         Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction,testBase,ansatzBase}
   *  \attention ret is assumed to be zero!
   */
  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
                const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunction, testBase, ansatzBase, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Codim0Interface< Traits, 2 >


/**
 *  \brief  Interface for local evaluations that depend on an intersection.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 */
template <class Traits, int numArguments>
class Codim1Interface
{
public:
  Codim1Interface() = delete;
};


/**
 *  \brief  Interface for binary codim 1 evaluations.
 */
template <class Traits>
class Codim1Interface<Traits, 2>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam L         Traits of the Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam T         Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A         Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction,testBase,ansatzBase}
   */
  template <class L, class T, class A, class D, int d, class R, int rL, int rCL, int rT, int rCT, int rA, int rCA>
  int order(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
            const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunction, testBase, ansatzBase));
    asImp().order(localFunction, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam LE        Traits of the entity Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam LN        Traits of the neighbor Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam TE        Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE        Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN        Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN        Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType Type of the codim 1 Intersection
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction*,testBase*,ansatzBase*}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction*,testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template <class L, class T, class A, class IntersectionType, class D, int d, class R, int rL, int rCL, int rT,
            int rCT, int rA, int rCA>
  void evaluate(const Dune::Stuff::LocalFunctionInterface<L, D, d, R, rL, rCL>& localFunction,
                const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase, const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate(localFunction, testBase, ansatzBase, intersection, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class class Codim1Interface< Traits, 2 >


/**
 *  \brief  Interface for quaternary codim 1 evaluations.
 */
template <class Traits>
class Codim1Interface<Traits, 4>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LE        Traits of the entity Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam LN        Traits of the neighbor Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam TE        Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE        Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN        Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN        Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction*,testBase*,ansatzBase*}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction*,testBase*,ansatzBase*}
   */
  template <class LE, class LN, class TE, class AE, class TN, class AN, class D, int d, class R, int rL, int rCL,
            int rT, int rCT, int rA, int rCA>
  int order(const Dune::Stuff::LocalFunctionInterface<LE, D, d, R, rL, rCL>& localFunctionEntity,
            const Dune::Stuff::LocalFunctionInterface<LN, D, d, R, rL, rCL>& localFunctionNeighbor,
            const BaseFunctionSetInterface<TE, D, d, R, rT, rCT>& testBaseEntity,
            const BaseFunctionSetInterface<AE, D, d, R, rA, rCA>& ansatzBaseEntity,
            const BaseFunctionSetInterface<TN, D, d, R, rT, rCT>& testBaseNeighbor,
            const BaseFunctionSetInterface<AN, D, d, R, rA, rCA>& ansatzBaseNeighbor) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctionEntity,
                                                 localFunctionNeighbor,
                                                 testBaseEntity,
                                                 ansatzBaseEntity,
                                                 testBaseNeighbor,
                                                 ansatzBaseNeighbor));
    asImp().order(localFunctionEntity,
                  localFunctionNeighbor,
                  testBaseEntity,
                  ansatzBaseEntity,
                  testBaseNeighbor,
                  ansatzBaseNeighbor);
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam LE        Traits of the entity Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam LN        Traits of the neighbor Dune::Stuff::LocalFunctionInterface implementation
   *  \tparam TE        Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE        Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN        Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN        Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType Type of the codim 1 Intersection
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction*,testBase*,ansatzBase*}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction*,testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template <class LE, class LN, class TE, class AE, class TN, class AN, class IntersectionType, class D, int d, class R,
            int rL, int rCL, int rT, int rCT, int rA, int rCA>
  void evaluate(const Dune::Stuff::LocalFunctionInterface<LE, D, d, R, rL, rCL>& localFunctionEntity,
                const Dune::Stuff::LocalFunctionInterface<LN, D, d, R, rL, rCL>& localFunctionNeighbor,
                const BaseFunctionSetInterface<TE, D, d, R, rT, rCT>& testBaseEntity,
                const BaseFunctionSetInterface<AE, D, d, R, rA, rCA>& ansatzBaseEntity,
                const BaseFunctionSetInterface<TN, D, d, R, rT, rCT>& testBaseNeighbor,
                const BaseFunctionSetInterface<AN, D, d, R, rA, rCA>& ansatzBaseNeighbor,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctionEntity,
                                                             localFunctionNeighbor,
                                                             testBaseEntity,
                                                             ansatzBaseEntity,
                                                             testBaseNeighbor,
                                                             ansatzBaseNeighbor,
                                                             intersection,
                                                             localPoint,
                                                             entityEntityRet,
                                                             neighborNeighborRet,
                                                             entityNeighborRet,
                                                             neighborEntityRet));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class class Codim1Interface< Traits, 4 >


} // namespace LocalEvaluation
} // namespace Discretizations
} // namespace Detailed
} // namespace Dune

#endif // DUNE_DETAILED_DISCRETIZATIONS_EVALUATION_INTERFACE_HH
