// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_EVALUATION_LAXWENDROFF_HH
#define DUNE_GDT_EVALUATION_LAXWENDROFF_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {
namespace LocalEvaluation {
namespace LaxWendroff {


// forwards
template< class LocalizableFunctionImp >
class Inner;

template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class Dirichlet;

template< class LocalizableFunctionImp >
class Absorbing;


namespace internal {


/**
 *  \brief  Traits for the Lax-Wendroff flux evaluation.
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
  typedef typename Dune::FieldVector< Dune::FieldMatrix< RangeFieldType, dimRange, dimRange >, dimDomain > FluxJacobianRangeType;
}; // class InnerTraits

/**
 *  \brief  Traits for the Lax-Wendroff flux evaluation at Dirichlet boundary intersections .
 */
template< class LocalizableFunctionImp, class BoundaryValueFunctionImp >
class DirichletTraits
    : public InnerTraits< LocalizableFunctionImp >
{
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
  typedef typename BoundaryValueFunctionType::DomainType  DomainType;
}; // class DirichletTraits

/**
 *  \brief  Traits for the Lax-Wendroff flux evaluation on absorbing boundary.
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
 *  \brief Lax-Wendroff flux evaluation for inner intersections and periodic boundary intersections.
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
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit Inner(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& ratio_dt_dx, const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , jacobian_(analytical_flux_.jacobian(RangeType(0)))
    , is_linear_(is_linear)
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
  void evaluate(const LocalfunctionTupleType& localFunctionsEntity,
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
    const auto local_center_neighbor = neighbor.geometry().local(neighbor.geometry().center());
    const std::vector< RangeType > u_j = ansatzBaseNeighbor.evaluate(local_center_neighbor);
    assert(u_i.size() == 1 && u_j.size() == 1);
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i[0]);
    const FluxRangeType f_u_j_temp = analytical_flux_.evaluate(u_j[0]);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_j;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
      f_u_j[ii] = f_u_j_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    size_t coord = 0;
    size_t num_zeros = 0;
    for (size_t ii = 0; ii < n_ij.size(); ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
        ++num_zeros;
      else
        DUNE_THROW(Dune::NotImplemented, "LaxWendroff flux is only implemented for axis parallel cube grids");
    }
    assert(num_zeros == dimRange - 1);
    FluxJacobianRangeType approx_jacobian = jacobian_;
    if (!is_linear_)
      reinitialize_jacobian(u_i[0], u_j[0], approx_jacobian);
    Dune::FieldVector< Dune::FieldVector< RangeFieldType, dimRange >, dimDomain > jacobian_multiplied;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      for (size_t kk = 0; kk < dimRange; ++kk) {
        for (size_t ll = 0; ll < dimRange; ++ll) {
          jacobian_multiplied[ii][kk] += approx_jacobian[ii][kk][ll]*(f_u_i[ll][ii] - f_u_j[ll][ii]);
        }
      }
    }
    RangeFieldType vol_intersection = 1;
    if (dimDomain != 1) {
      vol_intersection = intersection.geometry().volume();
    }
    const RangeFieldType ratio_dt_dx = std::get< 0 >(localFunctionsEntity)->evaluate(local_center_entity)[0];
    for (size_t kk = 0; kk < dimRange; ++kk)
      entityNeighborRet[kk][0] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - jacobian_multiplied[coord][kk]*ratio_dt_dx*0.5*n_ij[coord])*vol_intersection;
  } // void evaluate(...) const

private:
  void reinitialize_jacobian(const RangeType u_i,
                             const RangeType u_j,
                             FluxJacobianRangeType& jacobian) const
  {
    // clear jacobian
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian[ii] *= RangeFieldType(0);

    // calculate jacobian as jacobian(0.5*(u_i+u_j))
    RangeType u_mean = u_i;
    for (size_t ii = 0; ii < u_mean.size(); ++ii) {
      u_mean[ii] += u_j[ii];
      u_mean[ii] *= 0.5;
    }
    jacobian = analytical_flux_.jacobian(u_mean);
  } // void reinitialize_jacobian(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  FluxJacobianRangeType jacobian_;
  const bool is_linear_;
}; // class Inner

/**
 *  \brief  Lax-Wendroff flux evaluation for Dirichlet boundary intersections.
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
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::DomainType                               DomainType;
  typedef typename Traits::RangeType                                RangeType;
  static const unsigned int dimDomain = Traits::dimDomain;
  static const unsigned int dimRange = Traits::dimRange;

  // lambda = Delta t / Delta x
  explicit Dirichlet(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& ratio_dt_dx, const BoundaryValueFunctionType& boundary_values, const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
    , boundary_values_(boundary_values)
    , jacobian_(analytical_flux_.jacobian(RangeType(0)))
    , is_linear_(is_linear)
  {}

  LocalfunctionTupleType localFunctions(const EntityType& entity) const
  {
    return std::make_tuple(ratio_dt_dx_.local_function(entity), boundary_values_.local_function(entity));
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
  template< class IntersectionType, class R >
  void evaluate(const LocalfunctionTupleType& localFunctions,
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
    const auto local_center_intersection = entity.geometry().local(intersection.geometry().center());
    const auto& u_j = std::get< 1 >(localFunctions)->evaluate(local_center_intersection);
    assert(u_i.size() == 1);
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i[0]);
    const FluxRangeType f_u_j_temp = analytical_flux_.evaluate(u_j);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_j;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
      f_u_j[ii] = f_u_j_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    size_t coord = 0;
    size_t num_zeros = 0;
    for (size_t ii = 0; ii < n_ij.size(); ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
        ++num_zeros;
      else
        DUNE_THROW(Dune::NotImplemented, "LaxWendroff flux is only implemented for axis parallel cube grids");
    }
    assert(num_zeros == dimRange - 1);
    FluxJacobianRangeType approx_jacobian = jacobian_;
    if (!is_linear_)
      reinitialize_jacobian(u_i[0], u_j, approx_jacobian);
    Dune::FieldVector< Dune::FieldVector< RangeFieldType, dimRange >, dimDomain > jacobian_multiplied;
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      for (size_t kk = 0; kk < dimRange; ++kk) {
        for (size_t ll = 0; ll < dimRange; ++ll) {
          jacobian_multiplied[ii][kk] += approx_jacobian[ii][kk][ll]*(f_u_i[ll][ii] - f_u_j[ll][ii]);
        }
      }
    }
    RangeFieldType vol_intersection = 1;
    if (dimDomain != 1) {
      vol_intersection = intersection.geometry().volume();
    }
    const RangeFieldType ratio_dt_dx = std::get< 0 >(localFunctions)->evaluate(local_center_entity)[0];
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk][0] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - jacobian_multiplied[coord][kk]*ratio_dt_dx*0.5*n_ij[coord])*vol_intersection;
  } // void evaluate(...) const

private:
  void reinitialize_jacobian(const RangeType u_i,
                             const RangeType u_j,
                             FluxJacobianRangeType& jacobian) const
  {
    // clear jacobian
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian[ii] *= RangeFieldType(0);

    // calculate jacobian as jacobian(0.5*(u_i+u_j))
    RangeType u_mean = u_i;
    for (size_t ii = 0; ii < u_mean.size(); ++ii) {
      u_mean[ii] += u_j[ii];
      u_mean[ii] *= 0.5;
    }
    jacobian = analytical_flux_.jacobian(u_mean);
  } // void reinitialize_jacobian(...)

  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
  const BoundaryValueFunctionType& boundary_values_;
  FluxJacobianRangeType jacobian_;
  const bool is_linear_;
}; // class Dirichlet


/**
 *  \brief  Lax-Wendroff flux evaluation for absorbing boundary conditions on boundary intersections.
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

  explicit Absorbing(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& ratio_dt_dx, const bool = false)
    : analytical_flux_(analytical_flux)
    , ratio_dt_dx_(ratio_dt_dx)
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
    assert(u_i.size() == 1);
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i[0]);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(localPoint);
    RangeFieldType vol_intersection = 1;
    if (dimDomain != 1) {
      vol_intersection = intersection.geometry().volume();
    }
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk][0] = (f_u_i[kk] + f_u_i[kk])*n_ij*0.5*vol_intersection;
  } // void evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& ratio_dt_dx_;
}; // class Absorbing

} // namespace LaxWendroff
} // namespace Evaluation
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_EVALUATION_LAXWENDROFF_HH
