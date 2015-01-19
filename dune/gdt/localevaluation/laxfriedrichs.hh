// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Kirsten Weber, Tobias Leibner

#ifndef DUNE_GDT_EVALUATION_LAXFRIEDRICHS_HH
#define DUNE_GDT_EVALUATION_LAXFRIEDRICHS_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace LaxFriedrichs {


// forwards
template< class LocalizableFunctionImp >
class Inner;

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet;


namespace internal {


/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation.
 */
template< class LocalizableFunctionImp >
class InnerTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                                LocalizableFunctionType;
  typedef Inner< LocalizableFunctionType >                      derived_type;
  typedef typename LocalizableFunctionType::EntityType          EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType     DomainFieldType;
  typedef typename LocalizableFunctionType::RangeFieldType      RangeFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType   LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >    LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
  static const unsigned int dimRange = LocalizableFunctionType::dimRange;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  // copy of ProblemInterface::FluxType in dune/hdd/hyperbolic/problems/interfaces.hh
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity           FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType, RangeFieldType, dimRange, RangeFieldType, dimDomain > AnalyticalFluxType;
  typedef typename AnalyticalFluxType::RangeType                  FluxRangeType;
}; // class InnerTraits

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation at Dirichlet boundary intersections .
 */
template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class DirichletTraits
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, BoundaryValueFunctionImp >::value,
                "BoundaryValueFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                                LocalizableFunctionType;
  typedef BoundaryValueFunctionImp                              BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType BoundaryValueLocalfunctionType;
  typedef Dirichlet< LocalizableFunctionType, BoundaryValueFunctionType >                  derived_type;
  typedef typename LocalizableFunctionType::EntityType          EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType     DomainFieldType;
  typedef typename LocalizableFunctionType::RangeFieldType      RangeFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType   LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType >, std::shared_ptr< BoundaryValueLocalfunctionType > >    LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
  static const unsigned int dimRange = LocalizableFunctionType::dimRange;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  // copy of ProblemInterface::FluxType in dune/hdd/hyperbolic/problems/interfaces.hh
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity           FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType, RangeFieldType, dimRange, RangeFieldType, dimDomain > AnalyticalFluxType;
  typedef typename AnalyticalFluxType::RangeType                  FluxRangeType;
}; // class DirichletTraits
} // namespace internal


/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 */
template< class LocalizableFunctionImp >
class Inner
  : public LocalEvaluation::Codim1Interface< internal::InnerTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::InnerTraits< LocalizableFunctionImp >  Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  static const unsigned int dimDomain = Traits::dimDomain;

  // lambda = Delta t / Delta x
  explicit Inner(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& lambda)
    : analytical_flux_(analytical_flux)
    , lambda_(lambda)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(lambda_.local_function(entity));
  }

  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalfunctionTupleType& /*localFunctionsEntity*/,
               const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBaseNeighbor*/) const
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
                const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& /*entityEntityRet*/,
                Dune::DynamicMatrix< R >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< R >& entityNeighborRet,
                Dune::DynamicMatrix< R >& /*neighborEntityRet*/) const
  {
      const EntityType& entity = ansatzBaseEntity.entity();
      const EntityType& neighbor = ansatzBaseNeighbor.entity();

      const auto& u_i = ansatzBaseEntity.evaluate(entity.geometry().local(entity.geometry().center()));
      const auto& u_j = ansatzBaseNeighbor.evaluate(neighbor.geometry().local(neighbor.geometry().center()));
//      std::cout << "u_i 2 " << Dune::Stuff::Common::toString(u_i) <<" u_j " << Dune::Stuff::Common::toString(u_j) << std::endl;
      const FluxRangeType& f_u_i = analytical_flux_.evaluate(u_i[0]);
      const FluxRangeType& f_u_j = analytical_flux_.evaluate(u_j[0]);
//      std::cout << "f_u_i" << f_u_i << ", f_u_j" << f_u_j << std::endl;
      const auto local_center = entity.geometry().local(entity.geometry().center());
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      const auto lambda_ij = std::get< 0 >(localFunctionsEntity)->evaluate(local_center);
      std::cout << "ret:" << 1.0/2.0*(f_u_i*n_ij + f_u_j*n_ij) - 1.0/(2.0*lambda_ij[0])*(u_j[0] - u_i[0]) << std::endl;
      entityNeighborRet[0][0] =(f_u_i + f_u_j)*n_ij*0.5 - 1.0/(2.0*lambda_ij[0])*(u_j[0] - u_i[0]);
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& lambda_;
}; // class Inner

/**
 *  \brief  Lax-Friedrichs flux evaluation for Dirichlet boundary intersections.
 */
template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet
  : public LocalEvaluation::Codim1Interface
                          < internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >, 2 >
{
public:
  typedef internal::DirichletTraits< LocalizableFunctionImp, BoundaryValueFunctionImp >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  static const unsigned int dimDomain = Traits::dimDomain;

  // lambda = Delta t / Delta x
  explicit Dirichlet(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& lambda, const BoundaryValueFunctionType& boundary_values)
    : analytical_flux_(analytical_flux)
    , lambda_(lambda)
    , boundary_values_(boundary_values)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(lambda_.local_function(entity), boundary_values_.local_function(entity));
  }

  template< class R, int rT, int rCT, int rA, int rCA >
  size_t order(const LocalfunctionTupleType /*localFuncs*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBase*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& /*ansatzBase*/) const
  {
    DUNE_THROW(NotImplemented, "Not meant to be integrated");
  }

  /**
   *  \brief  Computes a binary codim 1 evaluation.
   *  \tparam IntersectionType    A model of Dune::Intersection< ... >
   *  \tparam R                   RangeFieldType
   *  \tparam r{T,A}              dimRange of the {testBase*,ansatzBase*}
   *  \tparam rC{T,A}             dimRangeRows of the {testBase*,ansatzBase*}
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R, int rT, int rCT, int rA, int rCA >
  void evaluate(const LocalfunctionTupleType& localFunctions,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rT, rCT >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, rA, rCA >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
      const EntityType& entity = ansatzBase.entity();
      const auto& u_i = ansatzBase.evaluate(entity.geometry().center());
      const auto local_intersection_center = entity.geometry().local(intersection.geometry().center());
      const auto& u_j = std::get< 1 >(localFunctions)->evaluate(local_intersection_center);
      //std::cout << "u_i 2 " << Dune::Stuff::Common::toString(u_i) <<" u_j " << Dune::Stuff::Common::toString(u_j) << std::endl;
      const FluxRangeType& f_u_i = analytical_flux_.evaluate(u_i[0]);
      const FluxRangeType& f_u_j = analytical_flux_.evaluate(u_j[0]);
      //std::cout << "f_u_i" << f_u_i << ", f_u_j" << f_u_j << std::endl;
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      const auto local_center = entity.geometry().local(entity.geometry().center());
      const auto lambda_ij = std::get< 0 >(localFunctions)->evaluate(local_center);
      //std::cout << "ret:" << 1.0/2.0*(f_u_i*n_ij + f_u_j*n_ij) - 1.0/(2.0*lambda_ij[0])*(u_j[0] - u_i[0]) << std::endl;
      ret[0][0] = (f_u_i + f_u_j)*n_ij*0.5 - 1.0/(2.0*lambda_ij[0])*(u_j[0] - u_i[0]);
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& lambda_;
  const BoundaryValueFunctionType& boundary_values_;
}; // class Dirichlet

} // namespace LaxFriedrichs
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_LAXFRIEDRICHS_HH
