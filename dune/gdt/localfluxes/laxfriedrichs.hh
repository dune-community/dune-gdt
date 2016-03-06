// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_LOCALFLUXES_LAXFRIEDRICHS_HH
#define DUNE_GDT_LOCALFLUXES_LAXFRIEDRICHS_HH

#include <tuple>
#include <memory>

# if HAVE_EIGEN
#   include <Eigen/Eigenvalues>
# endif

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/fmatrix.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/constant.hh>
#include <dune/stuff/la/container/eigen.hh>

#include "interfaces.hh"
#include "godunov.hh"

namespace Dune {
namespace GDT {


// forwards
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalCouplingFlux;

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalDirichletBoundaryFlux;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalAbsorbingBoundaryFlux;


namespace internal {


template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalCouplingFluxTraits
    : public GodunovNumericalCouplingFluxTraits< AnalyticalFluxImp, domainDim >
{
  static_assert(std::is_base_of< Dune::Stuff::IsLocalizableFunction, LocalizableFunctionImp >::value,
                "LocalizableFunctionImp has to be derived from Stuff::IsLocalizableFunction.");
public:
  typedef LocalizableFunctionImp                               LocalizableFunctionType;
  typedef typename LocalizableFunctionType::LocalfunctionType  LocalfunctionType;
  typedef std::tuple< std::shared_ptr< LocalfunctionType > >   LocalfunctionTupleType;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");

  typedef typename LocalizableFunctionType::DomainType              DomainType;
  typedef LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionType, domainDim >   derived_type;
}; // class LaxFriedrichsNumericalCouplingFluxTraits

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalDirichletBoundaryFluxTraits
    : public LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >
{
  typedef LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >    BaseType;
public:
  using BaseType::dimDomain;
  using typename BaseType::LocalfunctionType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType BoundaryValueLocalfunctionType;
  typedef LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >  derived_type;
  typedef std::tuple< std::shared_ptr< LocalfunctionType >,
                      std::shared_ptr< BoundaryValueLocalfunctionType > >  LocalfunctionTupleType;
}; // class LaxFriedrichsNumericalDirichletBoundaryFluxTraits

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalAbsorbingBoundaryFluxTraits
    : public LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >
{
public:
  typedef LaxFriedrichsNumericalAbsorbingBoundaryFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >  derived_type;
}; // class LaxFriedrichsNumericalAbsorbingBoundaryFluxTraits


} // namespace internal


/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class LaxFriedrichsNumericalCouplingFlux
  : public NumericalCouplingFluxInterface< internal::LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > >
{
public:
  typedef internal::LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  typedef typename LocalizableFunctionType::DomainType              DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                              const LocalizableFunctionType& dx,
                                              const double dt,
                                              const bool is_linear = false,
                                              const bool use_local = false,
                                              const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_local_(use_local)
    , entity_geometries_equal_(entity_geometries_equal)
  {
    geometry_evaluated_ = false;
  }

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
    RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));
    FluxRangeType f_u_i_plus_f_u_j = analytical_flux_.evaluate(u_i);
    f_u_i_plus_f_u_j += analytical_flux_.evaluate(u_j);
    auto n_ij = intersection.unitOuterNormal(x_intersection);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif //NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
#ifndef NDEBUG
        ++num_zeros;
#endif //NDEBUG
      }
      else
        DUNE_THROW(Dune::NotImplemented, "LaxFriedrichs flux is only implemented for axis parallel cube grids");
    }

    if (!use_local_) {
      const RangeFieldType dx = std::get< 0 >(local_functions_tuple_entity)->evaluate(intersection.geometryInInside().global(x_intersection))[0];
      *max_derivative_ = DomainType(dx/dt_);
    } else {
      if (!is_linear_ || !(*max_derivative_calculated_)) {
        *max_derivative_ = 0;
        const auto jacobian_u_i = analytical_flux_.jacobian(u_i);
        const auto jacobian_u_j = analytical_flux_.jacobian(u_j);
        std::vector< EigenMatrixType > jacobian_u_i_eigen;
        std::vector< EigenMatrixType > jacobian_u_j_eigen;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          jacobian_u_i_eigen.emplace_back(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_i[ii], 15)));
          jacobian_u_j_eigen.emplace_back(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_j[ii], 15)));
        }
      #if HAVE_EIGEN
          for (size_t ii = 0; ii < dimDomain; ++ii) {
            // create EigenSolver
            ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
                eigen_solver_u_i(jacobian_u_i_eigen[ii].backend());
            assert(eigen_solver_u_i.info() == ::Eigen::Success);
            ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
                eigen_solver_u_j(jacobian_u_j_eigen[ii].backend());
            assert(eigen_solver_u_j.info() == ::Eigen::Success);
            const auto eigenvalues_u_i = eigen_solver_u_i.eigenvalues(); // <- this should be an Eigen vector of std::complex
            assert(boost::numeric_cast< size_t >(eigenvalues_u_i.size()) == dimRange);
            const auto eigenvalues_u_j = eigen_solver_u_j.eigenvalues(); // <- this should be an Eigen vector of std::complex
            assert(boost::numeric_cast< size_t >(eigenvalues_u_j.size()) == dimRange);
            for (size_t jj = 0; jj < dimRange; ++jj) {
              // assert this is real
              assert(std::abs(eigenvalues_u_i[jj].imag()) < 1e-15);
              const RangeFieldType eigenvalue = eigenvalues_u_i[jj].real();
              if (std::abs(eigenvalue) > (*max_derivative_)[ii])
                (*max_derivative_)[ii] = std::abs(eigenvalue);
            }
            for (size_t jj = 0; jj < dimRange; ++jj) {
              // assert this is real
              assert(std::abs(eigenvalues_u_j[jj].imag()) < 1e-15);
              const RangeFieldType eigenvalue = eigenvalues_u_j[jj].real();
              if (std::abs(eigenvalue) > (*max_derivative_)[ii])
                (*max_derivative_)[ii] = std::abs(eigenvalue);
            }
          }
          if (is_linear_)
            *max_derivative_calculated_ = true;
      #else
          static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
      #endif
      }
    }
    if (!entity_geometries_equal_ || !(*geometry_evaluated_)) {
      *vol_intersection_ = intersection.geometry().volume();
      const auto& reference_element
        = Dune::ReferenceElements< DomainFieldType, dimDomain >::general(local_source_entity.entity().geometry().type());
      *num_neighbors_ = reference_element.size(1);
      *geometry_evaluated_ = true;
    }

    RangeType ret;
    //ret[kk] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - (u_j - u_i)[kk]*max_derivative_*1.0/num_neighbors_)*vol_intersection_
    //calculate (u_j - u_i)*max_derivative_/num_neighbors_*vol_intersection_
    u_j -= u_i;
    u_j *= (*max_derivative_)[coord]/(*num_neighbors_)*(*vol_intersection_);
    // scale n_ij by 0.5*vol_intersection_
    n_ij[coord] *= 0.5*vol_intersection_;
    // calculate flux
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = f_u_i_plus_f_u_j[kk][coord]*n_ij[coord] - u_j[kk];
    return ret;
  } // RangeType evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_local_;
  const bool entity_geometries_equal_;
  static typename DS::PerThreadValue< DomainType > max_derivative_;
  static typename DS::PerThreadValue< bool > max_derivative_calculated_;
  static typename DS::PerThreadValue< bool > geometry_evaluated_;
  mutable typename DS::PerThreadValue< RangeFieldType > vol_intersection_;
  mutable typename DS::PerThreadValue< int > num_neighbors_;
}; // class LaxFriedrichsNumericalCouplingFlux





template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< typename LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >::DomainType >
LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >::max_derivative_;

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >::max_derivative_calculated_(false);

template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >::geometry_evaluated_(false);

/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp >
class LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, 1>
  : public NumericalCouplingFluxInterface< internal::LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, 1 > >
{
public:
  typedef internal::LaxFriedrichsNumericalCouplingFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, 1 >     Traits;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                              const LocalizableFunctionType& dx,
                                              const double dt,
                                              const bool is_linear,
                                              const bool use_local = false,
                                              const bool /*entity_geometries_equal*/ = false)
    : analytical_flux_(analytical_flux)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_local_(use_local)
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
    RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));

    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    const RangeFieldType dx = std::get< 0 >(local_functions_tuple_entity)->evaluate(intersection.geometryInInside().global(x_intersection))[0];
    *max_derivative_ = dx/dt_;
    if (use_local_) {
      if (!is_linear_ || !(*max_derivative_calculated_)) {
        *max_derivative_ = 0;
        const auto jacobian_u_i = analytical_flux_.jacobian(u_i);
        const auto jacobian_u_j = analytical_flux_.jacobian(u_j);
        EigenMatrixType jacobian_u_i_eigen(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_i, 15)));
        EigenMatrixType jacobian_u_j_eigen(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_j, 15)));
      #if HAVE_EIGEN
        // create EigenSolver
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver_u_i(jacobian_u_i_eigen.backend());
        assert(eigen_solver_u_i.info() == ::Eigen::Success);
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver_u_j(jacobian_u_j_eigen.backend());
        assert(eigen_solver_u_j.info() == ::Eigen::Success);
        const auto eigenvalues_u_i = eigen_solver_u_i.eigenvalues(); // <- this should be an Eigen vector of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues_u_i.size()) == dimRange);
        const auto eigenvalues_u_j = eigen_solver_u_j.eigenvalues(); // <- this should be an Eigen vector of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues_u_j.size()) == dimRange);
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues_u_i[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues_u_i[jj].real();
          if (std::abs(eigenvalue) > *max_derivative_)
            *max_derivative_ = std::abs(eigenvalue);
        }
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues_u_j[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues_u_j[jj].real();
          if (std::abs(eigenvalue) > *max_derivative_)
            *max_derivative_ = std::abs(eigenvalue);
        }
        if (is_linear_)
          *max_derivative_calculated_ = true;
#else
        static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
      }
    }

    RangeType ret;
    // entityNeighborRet[0] = 0.5*((f(u_i) + f(u_j))*n_ij + max_derivative*(u_i - u_j)) where max_derivative = dx/dt if
    // we dont use the local LxF method. As the FieldVector does not provide an operator+, we have to split the expression.
    // calculate n_ij*(f(u_i) + f(u_j)) first
    ret = Dune::DynamicVector< RangeFieldType >(analytical_flux_.evaluate(u_i));
    ret += analytical_flux_.evaluate(u_j);
    if (n_ij < 0)
      ret *= n_ij;
    // add max_derivative*(u_i - u_j)
    u_i -= u_j;
    ret.axpy(*max_derivative_, u_i);
    // multiply by 0.5
    ret *= 0.5;
    return ret;
  } // void evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_local_;
  static typename DS::PerThreadValue< RangeFieldType > max_derivative_;
  static typename DS::PerThreadValue< bool >  max_derivative_calculated_;
}; // class LaxFriedrichsNumericalCouplingFlux< ... , 1 >

template< class AnalyticalFluxImp, class LocalizableFunctionImp >
typename DS::PerThreadValue< typename LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, 1 >::RangeFieldType >
LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, 1 >::max_derivative_(0);

template< class AnalyticalFluxImp, class LocalizableFunctionImp >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalCouplingFlux< AnalyticalFluxImp, LocalizableFunctionImp, 1 >::max_derivative_calculated_(false);



/**
*  \brief  Lax-Friedrichs flux evaluation for Dirichlet boundary intersections.
*/
template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain >
class LaxFriedrichsNumericalDirichletBoundaryFlux
    : public NumericalBoundaryFluxInterface< internal::LaxFriedrichsNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim > >
{
public:
  typedef internal::LaxFriedrichsNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::DomainType                               DomainType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsNumericalDirichletBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                       const BoundaryValueFunctionType& boundary_values,
                                                       const LocalizableFunctionType& dx,
                                                       const double dt,
                                                       const bool is_linear = false,
                                                       const bool use_local = false,
                                                       const bool entity_geometries_equal = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_local_(use_local)
    , entity_geometries_equal_(entity_geometries_equal)
  {
    geometry_evaluated_ = false;
  }

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
    FluxRangeType f_u_i_plus_f_u_j = analytical_flux_.evaluate(u_i);
    f_u_i_plus_f_u_j += analytical_flux_.evaluate(u_j);
    auto n_ij = intersection.unitOuterNormal(x_intersection);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif //NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
#ifndef NDEBUG
        ++num_zeros;
#endif //NDEBUG
      }
      else
        DUNE_THROW(Dune::NotImplemented, "Lax-Friedrichs flux is only implemented for axis parallel cube grids");
    }

    if (!use_local_) {
      const RangeFieldType dx = std::get< 0 >(local_functions_tuple)->evaluate(x_intersection_entity_coords)[0];
      *max_derivative_ = DomainType(dx/dt_);
    } else {
      if (!is_linear_ || !(*max_derivative_calculated_)) {
        *max_derivative_ = 0;
        const auto jacobian_u_i = analytical_flux_.jacobian(u_i);
        const auto jacobian_u_j = analytical_flux_.jacobian(u_j);
        std::vector< EigenMatrixType > jacobian_u_i_eigen;
        std::vector< EigenMatrixType > jacobian_u_j_eigen;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          jacobian_u_i_eigen.emplace_back(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_i[ii], 15)));
          jacobian_u_j_eigen.emplace_back(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_j[ii], 15)));
        }
      #if HAVE_EIGEN
          for (size_t ii = 0; ii < dimDomain; ++ii) {
            // create EigenSolver
            ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
                eigen_solver_u_i(jacobian_u_i_eigen[ii].backend());
            assert(eigen_solver_u_i.info() == ::Eigen::Success);
            ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
                eigen_solver_u_j(jacobian_u_j_eigen[ii].backend());
            assert(eigen_solver_u_j.info() == ::Eigen::Success);
            const auto eigenvalues_u_i = eigen_solver_u_i.eigenvalues(); // <- this should be an Eigen vector of std::complex
            assert(boost::numeric_cast< size_t >(eigenvalues_u_i.size()) == dimRange);
            const auto eigenvalues_u_j = eigen_solver_u_j.eigenvalues(); // <- this should be an Eigen vector of std::complex
            assert(boost::numeric_cast< size_t >(eigenvalues_u_j.size()) == dimRange);
            for (size_t jj = 0; jj < dimRange; ++jj) {
              // assert this is real
              assert(std::abs(eigenvalues_u_i[jj].imag()) < 1e-15);
              const RangeFieldType eigenvalue = eigenvalues_u_i[jj].real();
              if (std::abs(eigenvalue) > (*max_derivative_)[ii])
                (*max_derivative_)[ii] = std::abs(eigenvalue);
            }
            for (size_t jj = 0; jj < dimRange; ++jj) {
              // assert this is real
              assert(std::abs(eigenvalues_u_j[jj].imag()) < 1e-15);
              const RangeFieldType eigenvalue = eigenvalues_u_j[jj].real();
              if (std::abs(eigenvalue) > (*max_derivative_)[ii])
                (*max_derivative_)[ii] = std::abs(eigenvalue);
            }
          }
          if (is_linear_)
            *max_derivative_calculated_ = true;
      #else
          static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
      #endif
      }
    }
    if (!entity_geometries_equal_ || !(*geometry_evaluated_)) {
      *vol_intersection_ = intersection.geometry().volume();
      const auto& reference_element
        = Dune::ReferenceElements< DomainFieldType, dimDomain >::general(local_source_entity.entity().geometry().type());
      *num_neighbors_ = reference_element.size(1);
      *geometry_evaluated_ = true;
    }

    RangeType ret;
    //ret[kk] = ((f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - (u_j - u_i)[kk]*max_derivative_*1.0/num_neighbors_)*vol_intersection_
    //calculate (u_j - u_i)*max_derivative_/num_neighbors_*vol_intersection_
    u_j -= u_i;
    u_j *= (*max_derivative_)[coord]/(*num_neighbors_)*(*vol_intersection_);
    // scale n_ij by 0.5*vol_intersection_
    n_ij[coord] *= 0.5*vol_intersection_;
    // calculate flux
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = f_u_i_plus_f_u_j[kk][coord]*n_ij[coord] - u_j[kk];
    return ret;
  } // RangeType evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_local_;
  const bool entity_geometries_equal_;
  static typename DS::PerThreadValue< DomainType > max_derivative_;
  static typename DS::PerThreadValue< bool > max_derivative_calculated_;
  static typename DS::PerThreadValue< bool > geometry_evaluated_;
  mutable typename DS::PerThreadValue< RangeFieldType > vol_intersection_;
  mutable typename DS::PerThreadValue< int > num_neighbors_;
}; // class LaxFriedrichsNumericalDirichletBoundaryFlux

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< typename LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >::DomainType >
LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >::max_derivative_;

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >::max_derivative_calculated_(false);

template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, domainDim >::geometry_evaluated_(false);


/**
*  \brief  Lax-Friedrichs flux evaluation for LaxFriedrichsNumericalDirichletBoundaryFlux boundary intersections.
*/
template< class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
class LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 >
    : public NumericalBoundaryFluxInterface< internal::LaxFriedrichsNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 > >
{
public:
  typedef internal::LaxFriedrichsNumericalDirichletBoundaryFluxTraits< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 >  Traits;
  typedef typename Traits::BoundaryValueFunctionType                BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType                  LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType                   LocalfunctionTupleType;
  typedef typename Traits::EntityType                               EntityType;
  typedef typename Traits::DomainFieldType                          DomainFieldType;
  typedef typename Traits::RangeFieldType                           RangeFieldType;
  typedef typename Traits::AnalyticalFluxType                       AnalyticalFluxType;
  typedef typename Traits::FluxRangeType                            FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType                    FluxJacobianRangeType;
  typedef typename Traits::RangeType                                RangeType;
  typedef typename Traits::DomainType                               DomainType;
  typedef typename Traits::EigenMatrixType                          EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsNumericalDirichletBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                     const BoundaryValueFunctionType& boundary_values,
                     const LocalizableFunctionType& dx,
                     const double dt,
                     const bool is_linear = false,
                     const bool use_local = false,
                     const bool /*entity_geometries_equal*/ = false)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , dx_(dx)
    , dt_(dt)
    , is_linear_(is_linear)
    , use_local_(use_local)
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
    RangeType u_i = local_source_entity.evaluate(x_intersection_entity_coords);
    const auto u_j = std::get< 1 >(local_functions_tuple)->evaluate(x_intersection_entity_coords);

    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    const RangeFieldType dx = std::get< 0 >(local_functions_tuple)->evaluate(x_intersection_entity_coords)[0];
    *max_derivative_ = dx/dt_;
    if (use_local_) {
      if (!is_linear_ || !(*max_derivative_calculated_)) {
        *max_derivative_ = 0;
        const auto jacobian_u_i = analytical_flux_.jacobian(u_i);
        const auto jacobian_u_j = analytical_flux_.jacobian(u_j);
        EigenMatrixType jacobian_u_i_eigen(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_i, 15)));
        EigenMatrixType jacobian_u_j_eigen(DSC::from_string< EigenMatrixType >(DSC::to_string(jacobian_u_j, 15)));
      #if HAVE_EIGEN
        // create EigenSolver
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver_u_i(jacobian_u_i_eigen.backend());
        assert(eigen_solver_u_i.info() == ::Eigen::Success);
        ::Eigen::EigenSolver< typename Stuff::LA::EigenDenseMatrix< RangeFieldType >::BackendType >
            eigen_solver_u_j(jacobian_u_j_eigen.backend());
        assert(eigen_solver_u_j.info() == ::Eigen::Success);
        const auto eigenvalues_u_i = eigen_solver_u_i.eigenvalues(); // <- this should be an Eigen vector of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues_u_i.size()) == dimRange);
        const auto eigenvalues_u_j = eigen_solver_u_j.eigenvalues(); // <- this should be an Eigen vector of std::complex
        assert(boost::numeric_cast< size_t >(eigenvalues_u_j.size()) == dimRange);
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues_u_i[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues_u_i[jj].real();
          if (std::abs(eigenvalue) > *max_derivative_)
            *max_derivative_ = std::abs(eigenvalue);
        }
        for (size_t jj = 0; jj < dimRange; ++jj) {
          // assert this is real
          assert(std::abs(eigenvalues_u_j[jj].imag()) < 1e-15);
          const RangeFieldType eigenvalue = eigenvalues_u_j[jj].real();
          if (std::abs(eigenvalue) > *max_derivative_)
            *max_derivative_ = std::abs(eigenvalue);
        }
        if (is_linear_)
          *max_derivative_calculated_ = true;
#else
        static_assert(AlwaysFalse< FluxJacobianRangeType >::value, "You are missing eigen!");
#endif
      }
    }

    RangeType ret;
    // ret[0] = 0.5*((f(u_i) + f(u_j))*n_ij + max_derivative*(u_i - u_j)) where max_derivative = dx/dt if
    // we dont use the local LxF method. As the FieldVector does not provide an operator+, we have to split the expression.
    // calculate n_ij*(f(u_i) + f(u_j)) first
    ret = Dune::DynamicVector< RangeFieldType >(analytical_flux_.evaluate(u_i));
    ret += analytical_flux_.evaluate(u_j);
    if (n_ij < 0)
      ret *= n_ij;
    // add max_derivative*(u_i - u_j)
    u_i -= u_j;
    ret.axpy(*max_derivative_, u_i);
    // multiply by 0.5
    ret *= 0.5;
    return ret;
  } // RangeType evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  const LocalizableFunctionType& dx_;
  const double dt_;
  const bool is_linear_;
  const bool use_local_;
  static typename DS::PerThreadValue< RangeFieldType > max_derivative_;
  static typename DS::PerThreadValue< bool >  max_derivative_calculated_;
}; // class LaxFriedrichsNumericalDirichletBoundaryFlux< ... , 1 >

template < class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
typename DS::PerThreadValue< typename LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 >::RangeFieldType >
LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 >::max_derivative_(0);

template < class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp >
typename DS::PerThreadValue< bool >
LaxFriedrichsNumericalDirichletBoundaryFlux< AnalyticalFluxImp, BoundaryValueFunctionImp, LocalizableFunctionImp, 1 >::max_derivative_calculated_(false);



/**
 *  \brief  Lax-Friedrichs flux evaluation for absorbing boundary conditions on boundary intersections.
 */
template< class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim >
class LaxFriedrichsNumericalAbsorbingBoundaryFlux
  : public NumericalBoundaryFluxInterface< internal::LaxFriedrichsNumericalAbsorbingBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim > >
{
public:
  typedef internal::LaxFriedrichsNumericalAbsorbingBoundaryFluxTraits< AnalyticalFluxImp, LocalizableFunctionImp, domainDim >       Traits;
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

  explicit LaxFriedrichsNumericalAbsorbingBoundaryFlux(const AnalyticalFluxType& analytical_flux, const LocalizableFunctionType& /*dx*/, const double /*dt*/, const bool /*use_local*/)
    : analytical_flux_(analytical_flux)
  {}

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return std::make_tuple();
  }

  template< class IntersectionType >
  RangeType evaluate(const LocalfunctionTupleType& /*local_functions_tuple*/,
                      const Stuff::LocalfunctionInterface< EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1 >& local_source_entity,
                      const IntersectionType& intersection,
                      const Dune::FieldVector< DomainFieldType, dimDomain - 1 >& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const FluxRangeType f_u_i_temp = analytical_flux_.evaluate(u_i);
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
      ret[0][kk] = (f_u_i[kk] + f_u_i[kk])*n_ij*0.5*vol_intersection;
    return ret;
  } // void evaluate(...) const

private:
  const AnalyticalFluxType& analytical_flux_;
}; // class LaxFriedrichsNumericalAbsorbingBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFLUXES_LAXFRIEDRICHS_HH
