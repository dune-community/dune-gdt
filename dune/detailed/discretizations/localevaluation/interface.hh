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
 *
 *  A derived class has to provide a method
\code
template< class EntityType >
std::tuple< ... > localFunctions(const EntityType& entity) const
{
  ...
}
\endcode
 *  which returns the inducing local functions which are needed for one entity. Whatever this method returns will be
 *  given to order() and evaluate()
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
   *  \tparam LocalFunctionTuple Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T   Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template <class LocalFunctionTuple, class T, class D, int d, class R, int r, int rC>
  int order(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase));
    return asImp().order(localFunctions, testBase);
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam LocalFunctionTuple Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T   Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template <class LocalFunctionTuple, class T, class D, int d, class R, int r, int rC>
  void evaluate(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Codim0Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 0 evaluations.
 **/
template <class Traits>
class Codim0Interface<Traits, 2>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T       Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A} dimRangeRows of the {testBase,ansatzBase}
   */
  template <class LocalFunctionTuple, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  int order(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase, ansatzBase));
    return asImp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 0 evaluation.
   *  \tparam LocalFunctionTuple Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T         Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A         Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction,testBase,ansatzBase}
   *  \attention ret is assumed to be zero!
   */
  template <class LocalFunctionTuple, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase,
                const Dune::FieldVector<D, d>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, ansatzBase, localPoint, ret));
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
 *
 *  A derived class has to provide a method
\code
template< class EntityType >
std::tuple< ... > localFunctions(const EntityType& entity) const
{
  ...
}
\endcode
 *  which returns the inducing local functions which are needed for one entity. Whatever this method returns will be
 *  given to order() and evaluate();
 */
template <class Traits, int numArguments>
class Codim1Interface
{
public:
  Codim1Interface() = delete;
};


/**
 *  \brief  Interface for unary codim 1 evaluations.
 */
template <class Traits>
class Codim1Interface<Traits, 1>
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple  Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T                   Traits of the test BaseFunctionSetInterface implementation
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r                   dimRange of the testBase
   *  \tparam rC                  dimRangeRows of the testBase
   */
  template <class LocalFunctionTuple, class T, class D, int d, class R, int r, int rC>
  int order(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase));
    return asImp().order(localFunctions, testBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam LocalFunctionTuple  Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T                   Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r                   dimRange of the testBase
   *  \tparam rC                  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template <class LocalFunctionTuple, class T, class IntersectionType, class D, int d, class R, int r, int rC>
  void evaluate(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicVector<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, intersection, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class class Codim1Interface< Traits, 1 >


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
   *  \tparam LocalFunctionTuple  Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T                   Traits of the test BaseFunctionSetInterface implementation
   *  \tparam A                   Traits of the ansatz BaseFunctionSetInterface implementation
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase,ansatzBase}
   */
  template <class LocalFunctionTuple, class T, class A, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  int order(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
            const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase, ansatzBase));
    return asImp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam LocalFunctionTuple  Type of the localFunction container, that is returned by localFunctions(), i.e.
   * std::tuple< ... >
   *  \tparam T                   Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam A                   Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template <class LocalFunctionTuple, class T, class A, class IntersectionType, class D, int d, class R, int rT,
            int rCT, int rA, int rCA>
  void evaluate(const LocalFunctionTuple& localFunctions, const BaseFunctionSetInterface<T, D, d, R, rT, rCT>& testBase,
                const BaseFunctionSetInterface<A, D, d, R, rA, rCA>& ansatzBase, const IntersectionType& intersection,
                const Dune::FieldVector<D, d - 1>& localPoint, Dune::DynamicMatrix<R>& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(
        asImp().evaluate(localFunctions, testBase, ansatzBase, intersection, localPoint, ret));
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
   *  \tparam LocalFunctionTupleEn  Type of the entity localFunction container, that is returned by localFunctions(),
   * i.e. std::tuple< ... >
   *  \tparam LocalFunctionTupleNe  Type of the neighbor localFunction container, that is returned by localFunctions(),
   * i.e. std::tuple< ... >
   *  \tparam TE                    Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE                    Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN                    Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN                    Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam D                     DomainFieldType
   *  \tparam d                     dimDomain
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   */
  template <class LocalFunctionTupleEn, class LocalFunctionTupleNe, class TE, class AE, class TN, class AN, class D,
            int d, class R, int rT, int rCT, int rA, int rCA>
  int order(const LocalFunctionTupleEn& localFunctionsEntity, const LocalFunctionTupleNe& localFunctionsNeighbor,
            const BaseFunctionSetInterface<TE, D, d, R, rT, rCT>& testBaseEntity,
            const BaseFunctionSetInterface<AE, D, d, R, rA, rCA>& ansatzBaseEntity,
            const BaseFunctionSetInterface<TN, D, d, R, rT, rCT>& testBaseNeighbor,
            const BaseFunctionSetInterface<AN, D, d, R, rA, rCA>& ansatzBaseNeighbor) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctionsEntity,
                                                 localFunctionsNeighbor,
                                                 testBaseEntity,
                                                 ansatzBaseEntity,
                                                 testBaseNeighbor,
                                                 ansatzBaseNeighbor));
    return asImp().order(localFunctionsEntity,
                         localFunctionsNeighbor,
                         testBaseEntity,
                         ansatzBaseEntity,
                         testBaseNeighbor,
                         ansatzBaseNeighbor);
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam LocalFunctionTupleEn  Type of the entity localFunction container, that is returned by localFunctions(),
   * i.e. std::tuple< ... >
   *  \tparam LocalFunctionTupleNe  Type of the neighbor localFunction container, that is returned by localFunctions(),
   * i.e. std::tuple< ... >
   *  \tparam TE                    Traits of the entity test BaseFunctionSetInterface implementation
   *  \tparam AE                    Traits of the entity ansatz BaseFunctionSetInterface implementation
   *  \tparam TN                    Traits of the neighbor test BaseFunctionSetInterface implementation
   *  \tparam AN                    Traits of the neighbor ansatz BaseFunctionSetInterface implementation
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam D                     DomainFieldType
   *  \tparam d                     dimDomain
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template <class LocalFunctionTupleEn, class LocalFunctionTupleNe, class TE, class AE, class TN, class AN,
            class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA>
  void evaluate(const LocalFunctionTupleEn& localFunctionsEntity, const LocalFunctionTupleNe& localFunctionsNeighbor,
                const BaseFunctionSetInterface<TE, D, d, R, rT, rCT>& testBaseEntity,
                const BaseFunctionSetInterface<AE, D, d, R, rA, rCA>& ansatzBaseEntity,
                const BaseFunctionSetInterface<TN, D, d, R, rT, rCT>& testBaseNeighbor,
                const BaseFunctionSetInterface<AN, D, d, R, rA, rCA>& ansatzBaseNeighbor,
                const IntersectionType& intersection, const Dune::FieldVector<D, d - 1>& localPoint,
                Dune::DynamicMatrix<R>& entityEntityRet, Dune::DynamicMatrix<R>& neighborNeighborRet,
                Dune::DynamicMatrix<R>& entityNeighborRet, Dune::DynamicMatrix<R>& neighborEntityRet) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctionsEntity,
                                                             localFunctionsNeighbor,
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
