// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_EVALUATION_INTERFACE_HH
#define DUNE_GDT_EVALUATION_INTERFACE_HH

#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/bartonnackmanifcheck.hh>

#include <dune/stuff/functions/interfaces.hh>

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


/**
 *  \brief  Interface for local evaluations that depend on a codim 0 entity.
 *  \tparam numArguments  The number of local bases.
 *  \note   All evaluations have to be copyable!
 *
 *  A derived class has to provide a method
\code
template< class EntityType >
typename LocalfunctionTuple< EntityType >::Type localFunctions(const EntityType& entity) const
{
  ...
}
\endcode
 *  which returns the inducing local functions which are needed for one entity
 * and a class
\code
  template< class EntityType >
  class LocalfunctionTuple
  {
  public:
    typedef ... Type;
  };
\endcode
 * which defines the return type of this method. Whatever localFunctions() returns will be given to order() and
 * evaluate().
 */
template< class Traits, int numArguments >
class Codim0Interface
{
public:
  Codim0Interface() = delete;
};


/**
 *  \brief  Interface for unary codim 0 evaluations.
 */
template< class Traits >
class Codim0Interface< Traits, 1 >
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E   EntityType
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int r, int rC >
  size_t order(const LocalFunctionTuple& localFunctions,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase));
    return asImp().order(localFunctions, testBase);
  }

  /**
   *  \brief  Computes a unary codim 0 evaluation.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E   EntityType
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int r, int rC >
  void evaluate(const LocalFunctionTuple& localFunctions,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class Codim0Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 0 evaluations.
 **/
template< class Traits >
class Codim0Interface< Traits, 2 >
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E       EntityType
   *  \tparam D       DomainFieldType
   *  \tparam d       dimDomain
   *  \tparam R       RangeFieldType
   *  \tparam r{T,A}  dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A} dimRangeRows of the {testBase,ansatzBase}
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalFunctionTuple& localFunctions,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase, ansatzBase));
    return asImp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 0 evaluation.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E         EntityType
   *  \tparam D         DomainFieldType
   *  \tparam d         dimDomain
   *  \tparam R         RangeFieldType
   *  \tparam r{L,T,A}  dimRange of the {localFunction,testBase,ansatzBase}
   *  \tparam rC{L,T,A} dimRangeRows of the {localFunction,testBase,ansatzBase}
   *  \attention ret is assumed to be zero!
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalFunctionTuple& localFunctions,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                const Dune::FieldVector< D, d >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, ansatzBase, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
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
template< class Traits, int numArguments >
class Codim1Interface
{
public:
  Codim1Interface() = delete;
};


/**
 *  \brief  Interface for unary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 1 >
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E   EntityType
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange of the testBase
   *  \tparam rC  dimRangeRows of the testBase
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int r, int rC >
  size_t order(const LocalFunctionTuple& localFunctions,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase));
    return asImp().order(localFunctions, testBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam LocalFunctionTuple Type of the localFunction container that is returned by localFunctions()
   *  \tparam E                   EntityType
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r                   dimRange of the testBase
   *  \tparam rC                  dimRangeRows of the testBase
   *  \attention ret is assumed to be zero!
   */
  template< class LocalFunctionTuple, class E, class IntersectionType, class D, int d, class R, int r, int rC >
  void evaluate(const LocalFunctionTuple& localFunctions,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, r, rC >& testBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicVector< R >& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions, testBase, intersection, localPoint, ret));
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class class Codim1Interface< Traits, 1 >


/**
 *  \brief  Interface for binary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 2 >
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTuple  Type of the container that is returned by localFunctions()
   *  \tparam E                   EntityType
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase,ansatzBase}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase,ansatzBase}
   */
  template< class LocalFunctionTuple, class E, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalFunctionTuple& localFunctions,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctions, testBase, ansatzBase));
    return asImp().order(localFunctions, testBase, ansatzBase);
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam LocalFunctionTuple  Type of the container that is returned by localFunctions()
   *  \tparam E                   EntityType
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam D                   DomainFieldType
   *  \tparam d                   dimDomain
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template< class LocalFunctionTuple, class E, class IntersectionType,
            class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalFunctionTuple& localFunctions,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBase,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctions,
                                                             testBase, ansatzBase,
                                                             intersection,
                                                             localPoint,
                                                             ret));
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class class Codim1Interface< Traits, 2 >


/**
 *  \brief  Interface for quaternary codim 1 evaluations.
 */
template< class Traits >
class Codim1Interface< Traits, 4 >
{
public:
  typedef typename Traits::derived_type derived_type;

  /**
   *  \brief  Computes the needed integration order.
   *  \tparam LocalFunctionTupleEn  Type of the container that is returned by localFunctions(entity)
   *  \tparam LocalFunctionTupleNe  Type of the container that is returned by localFunctions(neighbour)
   *  \tparam E                     EntityType
   *  \tparam N                     NeighbourEntityType
   *  \tparam D                     DomainFieldType
   *  \tparam d                     dimDomain
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   */
  template< class LocalFunctionTupleEn, class LocalFunctionTupleNe, class E, class N,
            class D, int d, class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalFunctionTupleEn& localFunctionsEntity,
               const LocalFunctionTupleNe& localFunctionsNeighbor,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBaseEntity,
               const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBaseEntity,
               const Stuff::LocalfunctionSetInterface< N, D, d, R, rT, rCT >& testBaseNeighbor,
               const Stuff::LocalfunctionSetInterface< N, D, d, R, rA, rCA >& ansatzBaseNeighbor) const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().order(localFunctionsEntity, localFunctionsNeighbor,
                                                 testBaseEntity, ansatzBaseEntity,
                                                 testBaseNeighbor, ansatzBaseNeighbor));
    return asImp().order(localFunctionsEntity, localFunctionsNeighbor,
                         testBaseEntity, ansatzBaseEntity,
                         testBaseNeighbor, ansatzBaseNeighbor);
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam LocalFunctionTupleEn  Type of the container that is returned by localFunctions(entity)
   *  \tparam LocalFunctionTupleNe  Type of the container that is returned by localFunctions(neighbour)
   *  \tparam E                     EntityType
   *  \tparam N                     NeighbourEntityType
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam D                     DomainFieldType
   *  \tparam d                     dimDomain
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template< class LocalFunctionTupleEn, class LocalFunctionTupleNe, class E, class N,
            class IntersectionType, class D, int d, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalFunctionTupleEn& localFunctionsEntity,
                const LocalFunctionTupleNe& localFunctionsNeighbor,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rT, rCT >& testBaseEntity,
                const Stuff::LocalfunctionSetInterface< E, D, d, R, rA, rCA >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface< N, D, d, R, rT, rCT >& testBaseNeighbor,
                const Stuff::LocalfunctionSetInterface< N, D, d, R, rA, rCA >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< D, d - 1 >& localPoint,
                Dune::DynamicMatrix< R >& entityEntityRet,
                Dune::DynamicMatrix< R >& neighborNeighborRet,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& neighborEntityRet) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().evaluate(localFunctionsEntity,
                                                             localFunctionsNeighbor,
                                                             testBaseEntity, ansatzBaseEntity,
                                                             testBaseNeighbor, ansatzBaseNeighbor,
                                                             intersection,
                                                             localPoint,
                                                             entityEntityRet, neighborNeighborRet,
                                                             entityNeighborRet, neighborEntityRet));
  }

  const derived_type& asImp() const
  {
    return static_cast< const derived_type& >(*this);
  }
}; // class class Codim1Interface< Traits, 4 >


} // namespace LocalEvaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_INTERFACE_HH
