// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber, Tobias Leibner

#ifndef DUNE_GDT_EVALUATION_HYPERBOLIC_HH
#define DUNE_GDT_EVALUATION_HYPERBOLIC_HH

#include <tuple>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {


// forwards
template< class DiffusionFactorImp, class DiffusionTensorImp = void >
class LaxFriedrichsFlux;


namespace internal {


/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation.
 */
template< class LocalizableFunctionImp >
class LaxFriedrichsFluxTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                                LocalizableFunctionType;
  typedef LaxFriedrichsFlux< LocalizableFunctionType, void >    derived_type;
  typedef typename LocalizableFunctionType::EntityType          EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType     DomainFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType   LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >    LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
};
} // namespace internal


/**
 *  \brief  Computes an elliptic evaluation.
 */
template< class LocalizableFunctionImp >
class LaxFriedrichsFlux< LocalizableFunctionImp >
  : public LocalEvaluation::Codim1Interface< internal::LaxFriedrichsFluxTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::LaxFriedrichsFluxTraits< LocalizableFunctionImp, void >  Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  static const unsigned int dimDomain = Traits::dimDomain;

  explicit LaxFriedrichsFlux(const LocalizableFunctionType& analyticFlux)
    : analyticFlux_(analyticFlux)
    , dummyLocalizableFunction_(0)
  {}

  LocalfunctionTupleType dummyLocalFunction(const EntityType& entity) const
  {
    return std::make_tuple(dummyLocalizableFunction_.local_function(entity));
  }

  size_t order(const LocalfunctionTupleType localFunctionsEntity,
               const LocalfunctionTupleType localFunctionsNeighbor,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseEntity,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseEntity,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& testBaseNeighbor,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseNeighbor) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a quaternary codim 1 evaluation.
   *  \tparam IntersectionType      A model of Dune::Intersection< ... >
   *  \tparam R                     RangeFieldType
   *  \tparam r{T,A}                dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}               dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention entityEntityRet, entityEntityRet, entityEntityRet and neighborEntityRet are assumed to be zero!
   */
  template< class IntersectionType, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalfunctionTupleType& localFunctionsEntity,
                const LocalfunctionTupleType& localFunctionsNeighbor,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& /*localPoint*/,
                Dune::DynamicMatrix< R >& /*entityEntityRet*/,
                Dune::DynamicMatrix< R >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& /*neighborEntityRet*/) const
  {
    //TODO: implement
  }

  const LocalizableFunctionType& analyticFlux_;
  //dummy for the LocalOperator
  const Dune::Stuff::Functions::Constant dummyLocalizableFunction_;
}; // class LaxFriedrichsFlux


} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_HYPERBOLIC_HH
