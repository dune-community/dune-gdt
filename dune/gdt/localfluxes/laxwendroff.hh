// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCALFLUXES_LAXWENDROFF_HH
#define DUNE_GDT_LOCALFLUXES_LAXWENDROFF_HH

#include <tuple>
#include <memory>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>

#include "interfaces.hh"
#include "laxfriedrichs.hh"

namespace Dune {
namespace GDT {

// forwards
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxWendroffNumericalCouplingFlux;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
class LaxWendroffNumericalDirichletBoundaryFlux;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxWendroffNumericalAbsorbingBoundaryFlux;


namespace internal {


template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxWendroffNumericalCouplingFluxTraits
    : LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >
{
public:
  typedef LaxWendroffNumericalCouplingFlux< AnalyticalFluxImp,  LocalizableFunctionImp, domainDim >    derived_type;
}; // class LaxWendroffNumericalCouplingFluxTraits


template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim >
class LaxWendroffNumericalDirichletBoundaryFluxTraits
    : public LaxFriedrichsNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >
{
public:
  typedef LaxWendroffNumericalDirichletBoundaryFlux< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim >  derived_type;
}; // class LaxWendroffNumericalDirichletBoundaryFluxTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxWendroffNumericalAbsorbingBoundaryFluxTraits
    : public LaxWendroffNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >
{
public:
  typedef LaxWendroffNumericalAbsorbingBoundaryFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >  derived_type;
}; // class LaxWendroffNumericalAbsorbingBoundaryFluxTraits


} // namespace internal


/**
 *  \brief Lax-Wendroff flux evaluation for inner intersections and periodic boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class LaxWendroffNumericalCouplingFlux
  : public NumericalCouplingFluxInterface< internal::LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > >
{
public:
  typedef internal::LaxWendroffNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >           Traits;
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

  explicit LaxWendroffNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& dx, const double dt, const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , jacobian_(analytical_flux_.jacobian(RangeType(0)))
    , is_linear_(is_linear)
  {}

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity));
  }

  template< class IntersectionType >
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple_entity,
                      const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_neighbor,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i);
    const FluxRangeType f_u_j_temp = analytical_flux_.evaluate(u_j);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_j;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
      f_u_j[ii] = f_u_j_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
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
      reinitialize_jacobian(u_i, u_j, approx_jacobian);
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
    const RangeFieldType ratio_dt_dx = (std::get< 0 >(local_functions_tuple_entity)->evaluate(intersection.geometryInInside().global(x_intersection))[0])/dt_;
    RangeType ret;
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - jacobian_multiplied[coord][kk]*ratio_dt_dx*0.5*n_ij[coord])*vol_intersection;
    return ret;
  } // RangeType evaluate(...) const

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
  const LocalizableFunctionType& dx_;
  const double dt_;
  FluxJacobianRangeType jacobian_;
  const bool is_linear_;
}; // class LaxWendroffNumericalCouplingFlux

/**
 *  \brief  Lax-Wendroff flux evaluation for Dirichlet boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp, class BoundaryValueFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class LaxWendroffNumericalDirichletBoundaryFlux
    : public NumericalBoundaryFluxInterface< internal::LaxWendroffNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueFunctionImp, domainDim > >
{
public:
  typedef internal::LaxWendroffNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, BoundaryValueFunctionImp, 1 >  Traits;
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

  explicit LaxWendroffNumericalDirichletBoundaryFlux(const AnalyticalFluxType& analytical_flux, const BoundaryValueFunctionType& boundary_values, const LocalizableFunctionType& dx, const double dt, const bool is_linear = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , boundary_values_(boundary_values)
    , jacobian_(analytical_flux_.jacobian(RangeType(0)))
    , is_linear_(is_linear)
  {}

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(dx_.local_function(entity), boundary_values_.local_function(entity));
  }

  template< class IntersectionType >
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
    const RangeType u_i = local_source_entity.evaluate(x_intersection_entity_coords);
    auto u_j = std::get< 1 >(local_functions_tuple)->evaluate(x_intersection_entity_coords);
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i);
    const FluxRangeType f_u_j_temp = analytical_flux_.evaluate(u_j);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_j;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
      f_u_j[ii] = f_u_j_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
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
    const RangeFieldType ratio_dt_dx = (std::get< 0 >(local_functions_tuple)->evaluate(x_intersection_entity_coords))/dt_;
    RangeType ret;
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - jacobian_multiplied[coord][kk]*ratio_dt_dx*0.5*n_ij[coord])*vol_intersection;
    return ret;
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
  const LocalizableFunctionType& dx_;
  const double dt_;
  const BoundaryValueFunctionType& boundary_values_;
  FluxJacobianRangeType jacobian_;
  const bool is_linear_;
}; // class LaxWendroffNumericalDirichletBoundaryFlux


/**
 *  \brief  Lax-Wendroff flux evaluation for absorbing boundary conditions on boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxWendroffNumericalAbsorbingBoundaryFlux
  : public NumericalBoundaryFluxInterface< internal::LaxWendroffNumericalAbsorbingBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > >
{
public:
  typedef internal::LaxWendroffNumericalAbsorbingBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > Traits;
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

  explicit LaxWendroffNumericalAbsorbingBoundaryFlux(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& dx, const bool = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
  {}

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return std::make_tuple();
  }

  template< class IntersectionType >
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i);
    // copy to FieldMatrix (f_u_i_tmp is either a FieldVector or FieldMatrix)
    DSC::FieldMatrix< RangeFieldType, dimRange, dimDomain > f_u_i;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      f_u_i[ii] = f_u_i_temp[ii];
    }
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    RangeFieldType vol_intersection = 1;
    if (dimDomain != 1) {
      vol_intersection = intersection.geometry().volume();
    }
    RangeType ret;
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = (f_u_i[kk] + f_u_i[kk])*n_ij*0.5*vol_intersection;
    return ret;
  } // void evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
}; // class LaxWendroffNumericalAbsorbingBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFLUXES_LAXWENDROFF_HH
