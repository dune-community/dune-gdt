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

#include <dune/grid/yaspgrid.hh>

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

template< class LocalizableFunctionImp >
class Godunov;

template< class LocalizableFunctionImp >
class Absorbing;


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
  typedef LocalizableFunctionImp                                    LocalizableFunctionType;
  typedef Inner< LocalizableFunctionType >                          derived_type;
  typedef typename LocalizableFunctionType::EntityType              EntityType;
  typedef typename LocalizableFunctionType::DomainFieldType         DomainFieldType;
  typedef typename LocalizableFunctionType::RangeFieldType          RangeFieldType;
  typedef typename LocalizableFunctionType::LocalfunctionType       LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >        LocalfunctionTupleType;
  static const unsigned int dimDomain = LocalizableFunctionType::dimDomain;
  static const unsigned int dimRange = LocalizableFunctionType::dimRange;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  typedef typename Dune::YaspGrid< dimRange >::template Codim< 0 >::Entity              FluxSourceEntityType;
  typedef Dune::Stuff::GlobalFunctionInterface< FluxSourceEntityType,
                                                RangeFieldType, dimRange,
                                                RangeFieldType, dimRange, dimDomain >   AnalyticalFluxType;

  typedef typename AnalyticalFluxType::RangeType                                        FluxRangeType;
  typedef typename Dune::Stuff::LocalfunctionSetInterface< EntityType,
                                                           DomainFieldType, dimDomain,
                                                           RangeFieldType, dimRange, 1 >::RangeType  RangeType;
}; // class InnerTraits

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation at Dirichlet boundary intersections .
 */
template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class DirichletTraits
    : public InnerTraits< LocalizableFunctionImp >
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, BoundaryValueFunctionImp >::value,
                "BoundaryValueFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
  typedef InnerTraits< LocalizableFunctionImp >                            BaseType;
public:
  typedef LocalizableFunctionImp                                           LocalizableFunctionType;
  typedef BoundaryValueFunctionImp                                         BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType            BoundaryValueLocalfunctionType;
  typedef Dirichlet< LocalizableFunctionType, BoundaryValueFunctionType >  derived_type;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType >,
                      std::shared_ptr< BoundaryValueLocalfunctionType > >  LocalfunctionTupleType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RangeType;
}; // class DirichletTraits


template< class LocalizableFunctionImp >
class GodunovTraits
   : public InnerTraits< LocalizableFunctionImp >
{
  typedef InnerTraits< LocalizableFunctionImp > BaseType;
public:
  typedef LocalizableFunctionImp                LocalizableFunctionType;
  typedef Godunov< LocalizableFunctionType >    derived_type;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::LocalfunctionType;
  using typename BaseType::LocalfunctionTupleType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RangeType;
}; // class GodunovTraits

/**
 *  \brief  Traits for the Lax-Friedrichs flux evaluation on absorbing boundary.
 */
template< class LocalizableFunctionImp >
class AbsorbingTraits
   : public InnerTraits< LocalizableFunctionImp >
{
  typedef InnerTraits< LocalizableFunctionImp > BaseType;
public:
  typedef LocalizableFunctionImp                LocalizableFunctionType;
  typedef Absorbing< LocalizableFunctionType >  derived_type;
  using typename BaseType::EntityType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::RangeFieldType;
  using typename BaseType::LocalfunctionType;
  typedef typename std::tuple< >                LocalfunctionTupleType;
  using BaseType::dimDomain;
  using BaseType::dimRange;
  using typename BaseType::FluxSourceEntityType;
  using typename BaseType::AnalyticalFluxType;
  using typename BaseType::FluxRangeType;
  using typename BaseType::RangeType;
}; // class AbsorbingTraits


} // namespace internal


/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 */
template< class LocalizableFunctionImp >
class Inner
  : public LocalEvaluation::Codim1Interface< internal::InnerTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::InnerTraits< LocalizableFunctionImp >           Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::RangeType                                RangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit Inner(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& ratio_dt_dx)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(ratio_dt_dx_.local_function(entity));
  }

  size_t order(const LocalfunctionTupleType& /*localFunctionsEntity*/,
               const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*ansatzBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*ansatzBaseNeighbor*/) const
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
  template< class IntersectionType >
  void evaluate(const LocalfunctionTupleType& /*localFunctionsEntity*/,
                const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< RangeFieldType >& /*entityEntityRet*/,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< RangeFieldType >& entityNeighborRet,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborEntityRet*/) const
  {
    const EntityType& entity = ansatzBaseEntity.entity();
    const EntityType& neighbor = ansatzBaseNeighbor.entity();
    const auto local_center_entity = entity.geometry().local(entity.geometry().center());
    const std::vector< RangeType > u_i = ansatzBaseEntity.evaluate(local_center_entity);
    const std::vector< RangeType > u_j = ansatzBaseNeighbor.evaluate(neighbor.geometry().local(neighbor.geometry().center()));
    assert(u_i.size() == 1 && u_j.size() == 1);
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i[0]);
    const FluxRangeType f_u_j = analytical_flux_.evaluate(u_j[0]);
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    const auto jacobian_u_i = analytical_flux_.jacobian(u_i[0]);
    const auto jacobian_u_j = analytical_flux_.jacobian(u_j[0]);
    const auto max_derivative = jacobian_u_i.infinity_norm() > jacobian_u_j.infinity_norm()
                                          ? jacobian_u_i.infinity_norm()
                                          : jacobian_u_j.infinity_norm();
    for (size_t kk = 0; kk < dimRange; ++kk)
      entityNeighborRet[kk][0] = 0.5*((f_u_i[kk] + f_u_j[kk])*n_ij - (u_j[0] - u_i[0])[kk]*max_derivative);
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
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

  template< class R, unsigned long rT, unsigned long rCT, unsigned long rA, unsigned long rCA >
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
  template< class IntersectionType, class R, unsigned long rT, unsigned long rCT, unsigned long rA, unsigned long rCA >
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
      const auto local_center_entity = entity.geometry().local(entity.geometry().center());
      const auto& u_i = ansatzBase.evaluate(local_center_entity);
      const auto local_center_intersection = entity.geometry().local(intersection.geometry().center());
      const auto& u_j = std::get< 1 >(localFunctions)->evaluate(local_center_intersection);
      //std::cout << "u_i 2 " << Dune::Stuff::Common::toString(u_i) <<" u_j " << Dune::Stuff::Common::toString(u_j) << std::endl;
      const FluxRangeType& f_u_i = analytical_flux_.evaluate(u_i[0]);
      const FluxRangeType& f_u_j = analytical_flux_.evaluate(u_j[0]);
      //std::cout << "f_u_i" << f_u_i << ", f_u_j" << f_u_j << std::endl;
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      const auto lambda_ij = std::get< 0 >(localFunctions)->evaluate(local_center_entity);
      //std::cout << "ret:" << 1.0/2.0*(f_u_i*n_ij + f_u_j*n_ij) - 1.0/(2.0*lambda_ij[0])*(u_j[0] - u_i[0]) << std::endl;
      ret[0][0] = 0.5*((f_u_i + f_u_j)*n_ij - (u_j[0] - u_i[0])/lambda_ij[0]);
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& lambda_;
  const BoundaryValueFunctionType& boundary_values_;
}; // class Dirichlet



template< class LocalizableFunctionImp >
class Godunov
  : public LocalEvaluation::Codim1Interface< internal::GodunovTraits< LocalizableFunctionImp >, 4 >
{
public:
  typedef internal::GodunovTraits< LocalizableFunctionImp >           Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::AnalyticalFluxJacobianType             AnalyticalFluxJacobianType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                  FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  static const unsigned int dimDomain = Traits::dimDomain;

  explicit Godunov(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& ratio_dt_dx)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(ratio_dt_dx_.local_function(entity));
  }

  size_t order(const LocalfunctionTupleType& /*localFunctionsEntity*/,
               const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*ansatzBaseEntity*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseNeighbor*/,
               const Stuff::LocalfunctionSetInterface
                   < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*ansatzBaseNeighbor*/) const
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
  template< class IntersectionType >
  void evaluate(const LocalfunctionTupleType& /*localFunctionsEntity*/,
                const LocalfunctionTupleType& /*localFunctionsNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseEntity*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& ansatzBaseEntity,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& /*testBaseNeighbor*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, RangeFieldType, 1, 1 >& ansatzBaseNeighbor,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< RangeFieldType >& /*entityEntityRet*/,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborNeighborRet*/,
                Dune::DynamicMatrix< RangeFieldType >& entityNeighborRet,
                Dune::DynamicMatrix< RangeFieldType >& /*neighborEntityRet*/) const
  {
    const EntityType& entity = ansatzBaseEntity.entity();
    const EntityType& neighbor = ansatzBaseNeighbor.entity();
    const auto local_center_entity = entity.geometry().local(entity.geometry().center());
    const std::vector< RangeType > u_i = ansatzBaseEntity.evaluate(local_center_entity);
    const std::vector< RangeType > u_j = ansatzBaseNeighbor.evaluate(neighbor.geometry().local(neighbor.geometry().center()));
    //      std::cout << "u_i 2 " << Dune::Stuff::Common::toString(u_i) <<" u_j " << Dune::Stuff::Common::toString(u_j) << std::endl;
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i[0]);
    const FluxRangeType f_u_j = analytical_flux_.evaluate(u_j[0]);
    //      std::cout << "f_u_i" << f_u_i << ", f_u_j" << f_u_j << std::endl;
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    if (Dune::Stuff::Common::FloatCmp::eq(f_u_i, f_u_j)) {
      entityNeighborRet[0][0] = f_u_i*n_ij;
    } else {
      const RangeFieldType s = (f_u_i - f_u_j)/(u_i[0] - u_j[0]);
      if (s > 0 /*&& u_j[0] > 0*/) {
        entityNeighborRet[0][0] = 0.5*((f_u_j+f_u_i) - (f_u_j-f_u_i)*n_ij)*n_ij;
      } else if (s < 0 /*&& u_i[0] < 0*/) {
        entityNeighborRet[0][0] = 0.5*((f_u_j+f_u_i) - (f_u_i-f_u_j)*n_ij)*n_ij;
      } else {
        entityNeighborRet[0][0] = f_u_i*n_ij;
      }
    }
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
}; // class Godunov



/**
 *  \brief  Lax-Friedrichs flux evaluation for absorbing boundary conditions on boundary intersections.
 */
template< class LocalizableFunctionImp >
class Absorbing
  : public LocalEvaluation::Codim1Interface
                          < internal::AbsorbingTraits< LocalizableFunctionImp >, 2 >
{
public:
  typedef internal::AbsorbingTraits< LocalizableFunctionImp >       Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit Absorbing(const AnalyticalFluxType& analytical_flux)
    : analytical_flux_(analytical_flux)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& /*entity*/) const
  {
    return std::make_tuple();
  }

  template< class R, unsigned long rT, unsigned long rCT, unsigned long rA, unsigned long rCA >
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
   *  \attention ret is assumed to be zero!
   */
  template< class IntersectionType, class R >
  void evaluate(const LocalfunctionTupleType& /*localFunctions*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& /*testBase*/,
                const Stuff::LocalfunctionSetInterface
                    < EntityType, DomainFieldType, dimDomain, R, dimRange, 1 >& ansatzBase,
                const IntersectionType& intersection,
                const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& localPoint,
                Dune::DynamicMatrix< R >& ret) const
  {
      const EntityType& entity = ansatzBase.entity();
      const auto local_center_entity = entity.geometry().local(entity.geometry().center());
      const auto& u_i = ansatzBase.evaluate(local_center_entity);
      const FluxRangeType& f_u_i = analytical_flux_.evaluate(u_i[0]);
      const auto n_ij = intersection.unitOuterNormal(localPoint);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[kk][0] = 0.5*((f_u_i[kk] + f_u_i[kk])*n_ij);
  } // void evaluate(...) const

  const AnalyticalFluxType& analytical_flux_;
}; // class Absorbing

} // namespace LaxFriedrichs
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_LAXFRIEDRICHS_HH
