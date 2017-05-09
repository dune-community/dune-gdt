// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FLUXES_LAXFRIEDRICHS_HH
#define DUNE_GDT_LOCAL_FLUXES_LAXFRIEDRICHS_HH

#include <tuple>
#include <memory>

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <dune/common/typetraits.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/xt/functions/interfaces.hh>

#include "interfaces.hh"
#include "godunov.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalDirichletNumericalBoundaryFlux;

template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux;

#if HAVE_EIGEN

template <class XtLaMatrixImp, size_t rows, size_t cols>
struct FieldMatrixToLaDenseMatrixConverter
{
  typedef XtLaMatrixImp LaMatrixType;
  typedef Dune::FieldMatrix<typename LaMatrixType::ScalarType, rows, cols> FieldMatrixType;

  static LaMatrixType convert(const FieldMatrixType& in)
  {
    LaMatrixType out(rows, cols);
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        out.set_entry(ii, jj, in[ii][jj]);
    return out;
  }

  static FieldMatrixType convert_back(const LaMatrixType& in)
  {
    FieldMatrixType out;
    for (size_t ii = 0; ii < rows; ++ii)
      for (size_t jj = 0; jj < cols; ++jj)
        out[ii][jj] = in.get_entry(ii, jj);
    return out;
  }
};


namespace internal {


template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalNumericalCouplingFluxTraits
    : public GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, domainDim>
{
  static_assert(Dune::XT::Functions::is_localizable_function<LocalizableFunctionImp>::value,
                "LocalizableFunctionImp has to be derived from XT::Functions::is_localizable_function.");

public:
  typedef LocalizableFunctionImp LocalizableFunctionType;
  typedef typename LocalizableFunctionType::LocalfunctionType LocalfunctionType;
  typedef std::tuple<int, std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
  static_assert(LocalizableFunctionType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");

  typedef typename LocalizableFunctionType::DomainType DomainType;
  typedef LaxFriedrichsLocalNumericalCouplingFlux<AnalyticalFluxImp, LocalizableFunctionType, domainDim> derived_type;
}; // class LaxFriedrichsLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueFunctionImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits
    : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, domainDim>
{
  typedef LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, domainDim> BaseType;

public:
  using typename BaseType::LocalfunctionType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType BoundaryValueLocalfunctionType;
  typedef LaxFriedrichsLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                           BoundaryValueFunctionImp,
                                                           LocalizableFunctionImp,
                                                           domainDim>
      derived_type;
  typedef std::tuple<size_t, std::shared_ptr<LocalfunctionType>, std::shared_ptr<BoundaryValueLocalfunctionType>>
      LocalfunctionTupleType;
}; // class LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits

template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalAbsorbingNumericalBoundaryFluxTraits
    : public LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, domainDim>
{
public:
  typedef LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux<AnalyticalFluxImp, LocalizableFunctionImp, domainDim>
      derived_type;
}; // class LaxFriedrichsLocalAbsorbingNumericalBoundaryFluxTraits

template <class Traits>
class LaxFriedrichsFluxImplementation
{
public:
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsFluxImplementation(const AnalyticalFluxType& analytical_flux,
                                           const XT::Common::Parameter param,
                                           const bool use_local = false,
                                           const bool is_linear = false,
                                           const DomainType lambda = DomainType(0))
    : analytical_flux_(analytical_flux)
    , dt_(param.get("dt")[0])
    , t_(param.get("t")[0])
    , use_local_(use_local)
    , is_linear_(is_linear)
    , lambda_(lambda)
    , lambda_provided_(XT::Common::FloatCmp::ne(lambda_, DomainType(0)))
  {
    if (lambda_provided_ && use_local_)
      std::cerr << "WARNING: Parameter lambda in Lax-Friedrichs flux set but will have no effect because local "
                   "Lax-Friedrichs flux is requested."
                << std::endl;

    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
  }

  ~LaxFriedrichsFluxImplementation()
  {
    is_instantiated_ = false;
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                     const EntityType& inside_entity,
                     const EntityType& outside_entity,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords,
                     const DomainType& x_in_inside_coords,
                     const DomainType& x_in_outside_coords,
                     const RangeType& u_i,
                     const RangeType& u_j) const
  {
    auto f_u_i_plus_f_u_j = XT::Functions::internal::RangeTypeConverter<dimRange, dimDomain>::convert(
        analytical_flux_.evaluate(u_i, inside_entity, x_in_inside_coords, t_));
    f_u_i_plus_f_u_j += XT::Functions::internal::RangeTypeConverter<dimRange, dimDomain>::convert(
        analytical_flux_.evaluate(u_j, outside_entity, x_in_outside_coords, t_));
    auto n_ij = intersection.unitOuterNormal(x_in_intersection_coords);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif // NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (Dune::XT::Common::FloatCmp::eq(n_ij[ii], RangeFieldType(1))
          || Dune::XT::Common::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
      else if (Dune::XT::Common::FloatCmp::eq(n_ij[ii], RangeFieldType(0))) {
#ifndef NDEBUG
        ++num_zeros;
#endif // NDEBUG
      } else
        DUNE_THROW(Dune::NotImplemented, "LaxFriedrichs flux is only implemented for axis parallel cube grids");
    }

    if (use_local_) {
      if (!is_linear_ || !max_derivative_calculated_) {
        DomainType max_derivative(0);
        const auto jacobian_u_i =
            XT::Functions::internal::JacobianRangeTypeConverter<dimRange, dimRange, dimDomain>::convert(
                analytical_flux_.jacobian(u_i, inside_entity, x_in_inside_coords, t_));
        const auto jacobian_u_j =
            XT::Functions::internal::JacobianRangeTypeConverter<dimRange, dimRange, dimDomain>::convert(
                analytical_flux_.jacobian(u_j, outside_entity, x_in_outside_coords, t_));
        std::vector<EigenMatrixType> jacobian_u_i_eigen;
        std::vector<EigenMatrixType> jacobian_u_j_eigen;
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          jacobian_u_i_eigen.emplace_back(
              FieldMatrixToLaDenseMatrixConverter<EigenMatrixType, dimRange, dimRange>::convert(jacobian_u_i[ii]));
          jacobian_u_j_eigen.emplace_back(
              FieldMatrixToLaDenseMatrixConverter<EigenMatrixType, dimRange, dimRange>::convert(jacobian_u_j[ii]));
        }
        for (size_t ii = 0; ii < dimDomain; ++ii) {
          // create EigenSolver
          ::Eigen::EigenSolver<typename XT::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver_u_i(
              jacobian_u_i_eigen[ii].backend());
          assert(eigen_solver_u_i.info() == ::Eigen::Success);
          ::Eigen::EigenSolver<typename XT::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver_u_j(
              jacobian_u_j_eigen[ii].backend());
          assert(eigen_solver_u_j.info() == ::Eigen::Success);
          const auto eigenvalues_u_i =
              eigen_solver_u_i.eigenvalues(); // <- this should be an Eigen vector of std::complex
          assert(boost::numeric_cast<size_t>(eigenvalues_u_i.size()) == dimRange);
          const auto eigenvalues_u_j =
              eigen_solver_u_j.eigenvalues(); // <- this should be an Eigen vector of std::complex
          assert(boost::numeric_cast<size_t>(eigenvalues_u_j.size()) == dimRange);
          for (size_t jj = 0; jj < dimRange; ++jj) {
            // assert this is real
            assert(std::abs(eigenvalues_u_i[jj].imag()) < 1e-15);
            const RangeFieldType eigenvalue = eigenvalues_u_i[jj].real();
            if (std::abs(eigenvalue) > max_derivative[ii])
              max_derivative[ii] = std::abs(eigenvalue);
          }
          for (size_t jj = 0; jj < dimRange; ++jj) {
            // assert this is real
            assert(std::abs(eigenvalues_u_j[jj].imag()) < 1e-15);
            const RangeFieldType eigenvalue = eigenvalues_u_j[jj].real();
            if (std::abs(eigenvalue) > max_derivative[ii])
              max_derivative[ii] = std::abs(eigenvalue);
          }
        }
        if (is_linear_)
          max_derivative_calculated_ = true;
        for (size_t ii = 0; ii < dimDomain; ++ii)
          lambda_ij_[ii] = 1. / max_derivative[ii];
      }
    } else if (lambda_provided_) {
      lambda_ij_ = DomainType(lambda_);
    } else {
      const RangeFieldType dx = std::get<1>(local_functions_tuple)->evaluate(x_in_inside_coords)[0];
      lambda_ij_ = DomainType(dt_ / dx);
    } // if (use_local)

    // calculate flux evaluation as
    // ret[kk] = (f_u_i[kk] + f_u_j[kk])*n_ij*0.5 - (u_j - u_i)[kk]*1.0/(num_neighbors_*lambda_ij)
    RangeType ret;
    const size_t num_neighbors = std::get<0>(local_functions_tuple);
    auto second_part = u_j;
    second_part -= u_i;
    second_part /= lambda_ij_[coord] * num_neighbors;
    n_ij[coord] *= 0.5;
    for (size_t kk = 0; kk < dimRange; ++kk)
      ret[kk] = f_u_i_plus_f_u_j[kk][coord] * n_ij[coord] - second_part[kk];
    return ret;
  }

private:
  const AnalyticalFluxType& analytical_flux_;
  const double dt_;
  const double t_;
  const bool use_local_;
  const bool is_linear_;
  const DomainType lambda_;
  const bool lambda_provided_;
  static thread_local DomainType lambda_ij_;
  static thread_local bool max_derivative_calculated_;
  static thread_local int num_neighbors_;
  static bool is_instantiated_;
}; // class LaxFriedrichsFluxImplementation<...>

template <class Traits>
thread_local
    typename LaxFriedrichsFluxImplementation<Traits>::DomainType LaxFriedrichsFluxImplementation<Traits>::lambda_ij_;

template <class Traits>
thread_local bool LaxFriedrichsFluxImplementation<Traits>::max_derivative_calculated_(false);

template <class Traits>
bool LaxFriedrichsFluxImplementation<Traits>::is_instantiated_(false);


} // namespace internal


/**
 *  \brief  Lax-Friedrichs flux evaluation for inner intersections and periodic boundary intersections.
 *
 *  The Lax-Friedrichs flux is an approximation to the integral
 *  \int_{S_{ij}} \mathbf{F}(\mathbf{u}) \cdot \mathbf{n}_{ij},
 *  where S_{ij} is the intersection between the entities i and j, \mathbf{F}(\mathbf{u}) is the analytical flux
 *  (evaluated at \mathbf{u}) and \mathbf{n}_{ij} is the unit outer normal of S_{ij}.
 *  The Lax-Friedrichs flux takes the form
 *  \mathbf{g}_{ij}^{LF}(\mathbf{u}_i, \mathbf{u}_j)
 *  = \int_{S_{ij}} \frac{1}{2}(\mathbf{F}(\mathbf{u}_i) + \mathbf{F}(\mathbf{u}_j) \cdot \mathbf{n}_{ij}
 *  - \frac{1}{\alpha_i \lambda_{ij}} (\mathbf{u}_j - \mathbf{u}_i),
 *  where \alpha_i is the number of neighbors (i.e. intersections) of the entity i and lambda_{ij} is a local
 *  constant fulfilling
 *  \lambda_{ij} \sup_{\mathbf{u}} (\mathbf{F}(\mathbf{u} \cdot \mathbf{n}_{ij})^\prime \leq 1.
 *  The integration is done numerically and implemented in the LocalCouplingFvOperator. This class implements
 *  the evaluation of the integrand. As we are restricting ourselves to axis-parallel cubic grids, only one component of
 *  \mathbf{n}_{ij} is non-zero, denote this component by k. Then the Lax-Friedrichs flux evaluation reduces to
 *  \frac{1}{2}(\mathbf{f}^k(\mathbf{u}_i) + \mathbf{f}^k(\mathbf{u}_j) n_{ij,k}
 *  - \frac{1}{\alpha_i \lambda_{ij}} (\mathbf{u}_j - \mathbf{u}_i),
 *  where \mathbf{f}^k is the k-th column of the analytical flux.
 *  For the classical Lax-Friedrichs flux, \lambda_{ij} is chosen as dt/dx_i, where dt is the current time
 *  step length and dx_i is the width of entity i. This fulfills the equation above as long as the CFL condition
 *  is fulfilled.
 *  The local Lax-Friedrichs flux can be chosen by setting \param use_local to true, here \lambda_{ij} is chosen
 *  as the inverse of the maximal eigenvalue of \mathbf{f}^k(\mathbf{u}_i) and \mathbf{f}^k(\mathbf{u}_j). In this
 *  case, you should also specify whether your analytical flux is linear by setting \param is_linear, which avoids
 *  recalculating the eigenvalues on every intersection in the linear case.
 *  You can also provide a user-defined \param lambda that is used as \lambda_{ij} on all intersections. You need to set
 *  use_local to false, otherwise lambda will not be used.
 */
template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim = LocalizableFunctionImp::dimDomain>
class LaxFriedrichsLocalNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::
                                                     LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                   LocalizableFunctionImp,
                                                                                                   domainDim>>
{
public:
  typedef internal::LaxFriedrichsLocalNumericalCouplingFluxTraits<AnalyticalFluxImp, LocalizableFunctionImp, domainDim>
      Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  typedef typename LocalizableFunctionType::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                                   const LocalizableFunctionType& dx,
                                                   const XT::Common::Parameter& param,
                                                   const bool use_local = false,
                                                   const bool is_linear = false,
                                                   const DomainType lambda = DomainType(0))
    : dx_(dx)
    , implementation_(analytical_flux, param, use_local, is_linear, lambda)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(entity.subEntities(1), dx_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple_entity,
      const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_neighbor,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    const auto& inside_entity = intersection.inside();
    const auto& outside_entity = intersection.outside();
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const auto x_in_outside_coords = intersection.geometryInOutside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = local_source_neighbor.evaluate(x_in_outside_coords);
    return implementation_.evaluate(local_functions_tuple_entity,
                                    inside_entity,
                                    outside_entity,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_outside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const LocalizableFunctionType& dx_;
  const internal::LaxFriedrichsFluxImplementation<Traits> implementation_;
}; // class LaxFriedrichsLocalNumericalCouplingFlux

/**
*  \brief  Lax-Friedrichs flux evaluation for Dirichlet boundary intersections.
*  \see    LaxFriedrichsLocalNumericalCouplingFlux
*/
template <class AnalyticalFluxImp,
          class BoundaryValueFunctionImp,
          class LocalizableFunctionImp,
          size_t domainDim = LocalizableFunctionImp::dimDomain>
class LaxFriedrichsLocalDirichletNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                            BoundaryValueFunctionImp,
                                                                                                            LocalizableFunctionImp,
                                                                                                            domainDim>>
{
public:
  typedef internal::LaxFriedrichsLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                           BoundaryValueFunctionImp,
                                                                           LocalizableFunctionImp,
                                                                           domainDim>
      Traits;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsLocalDirichletNumericalBoundaryFlux(
      const AnalyticalFluxType& analytical_flux,
      const std::shared_ptr<BoundaryValueFunctionType>& boundary_values,
      const LocalizableFunctionType& dx,
      const XT::Common::Parameter param,
      const bool use_local = false,
      const bool is_linear = false,
      const DomainType lambda = DomainType(0))
    : boundary_values_(boundary_values)
    , dx_(dx)
    , implementation_(analytical_flux, param, use_local, is_linear, lambda)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(entity.subEntities(1), dx_.local_function(entity), boundary_values_->local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const

  {
    const auto& entity = intersection.inside();
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const auto x_in_outside_coords = DomainType(200);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = std::get<2>(local_functions_tuple)->evaluate(x_in_inside_coords);
    return implementation_.evaluate(local_functions_tuple,
                                    entity,
                                    entity,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_outside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

private:
  const std::shared_ptr<BoundaryValueFunctionType>& boundary_values_;
  const LocalizableFunctionType& dx_;
  const internal::LaxFriedrichsFluxImplementation<Traits> implementation_;
}; // class LaxFriedrichsLocalDirichletNumericalBoundaryFlux

/**
 *  \brief  Lax-Friedrichs flux evaluation for absorbing boundary conditions on boundary intersections.
 *  \see    LaxFriedrichsLocalNumericalCouplingFlux
 */
template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     LaxFriedrichsLocalAbsorbingNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                            LocalizableFunctionImp,
                                                                                                            domainDim>>
{
public:
  typedef internal::LaxFriedrichsLocalAbsorbingNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                           LocalizableFunctionImp,
                                                                           domainDim>
      Traits;
  typedef typename Traits::LocalizableFunctionType LocalizableFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                            const LocalizableFunctionType& dx,
                                                            const XT::Common::Parameter param,
                                                            const bool use_local = false,
                                                            const bool is_linear = false,
                                                            const DomainType lambda = DomainType(0))
    : dx_(dx)
    , implementation_(analytical_flux, param, use_local, is_linear, lambda)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(entity.subEntities(1), dx_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple_entity,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const
  {
    const auto& entity = intersection.inside();
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    return implementation_.evaluate(local_functions_tuple_entity,
                                    entity,
                                    entity,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_inside_coords,
                                    u_i,
                                    u_i);
  } // void evaluate(...) const

private:
  const LocalizableFunctionType& dx_;
  const internal::LaxFriedrichsFluxImplementation<Traits>& implementation_;
}; // class LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux


#else // HAVE_EIGEN

template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalNumericalCouplingFlux
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionType, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalDirichletNumericalBoundaryFlux
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class LocalizableFunctionImp, size_t domainDim>
class LaxFriedrichsLocalAbsorbingNumericalBoundaryFlux
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_LAXFRIEDRICHS_HH
