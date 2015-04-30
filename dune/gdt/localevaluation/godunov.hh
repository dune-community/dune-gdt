// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_EVALUATION_GODUNOV_HH
#define DUNE_GDT_EVALUATION_GODUNOV_HH

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
namespace Godunov {

template< class LocalizableFunctionImp >
class Inner;

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


} // namespace internal


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


} // namespace Godunov
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_GODUNOV_HH
