// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
#define DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH

#include <tuple>
#include <memory>

#if HAVE_EIGEN
#include <Eigen/Eigenvalues>
#endif

#include <boost/numeric/conversion/cast.hpp>

#include <dune/common/dynmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/stuff/common/string.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/functions/affine.hh>
#include <dune/stuff/la/container/eigen.hh>
#include <dune/stuff/common/parallel/threadstorage.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, size_t domainDim>
class LocalGodunovNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueFunctionType, size_t domainDim>
class LocalGodunovNumericalBoundaryFlux;

// TODO: remove Eigen-dependency and use generic eigenvalue solver
#if HAVE_EIGEN


namespace internal {


template <class AnalyticalFluxImp, size_t domainDim = AnalyticalFluxImp::dimDomain>
class LocalGodunovNumericalCouplingFluxTraits
{
  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef LocalGodunovNumericalCouplingFlux<AnalyticalFluxType, domainDim> derived_type;
  typedef typename AnalyticalFluxType::EntityType EntityType;
  typedef typename AnalyticalFluxType::DomainFieldType DomainFieldType;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename AnalyticalFluxType::RangeType RangeType;
  typedef typename AnalyticalFluxType::FluxRangeType FluxRangeType;
  typedef typename AnalyticalFluxType::FluxJacobianRangeType FluxJacobianRangeType;
  typedef std::tuple<double> LocalfunctionTupleType;
  static const size_t dimDomain = domainDim;
  static const size_t dimRange  = AnalyticalFluxType::dimRange;
  static_assert(AnalyticalFluxType::dimRangeCols == 1, "Not implemented for dimRangeCols > 1!");
  typedef typename Dune::Stuff::LA::EigenDenseMatrix<RangeFieldType> EigenMatrixType;
}; // class LocalGodunovNumericalCouplingFluxTraits

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp,
          size_t domainDim = AnalyticalBoundaryFluxImp::dimDomain>
class LocalGodunovNumericalBoundaryFluxTraits
    : public LocalGodunovNumericalCouplingFluxTraits<AnalyticalBoundaryFluxImp, domainDim>
{
  typedef LocalGodunovNumericalCouplingFluxTraits<AnalyticalBoundaryFluxImp, domainDim> BaseType;

public:
  typedef AnalyticalBoundaryFluxImp AnalyticalFluxType;
  typedef BoundaryValueFunctionImp BoundaryValueFunctionType;
  typedef typename BoundaryValueFunctionType::LocalfunctionType LocalfunctionType;
  typedef LocalGodunovNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueFunctionType, domainDim> derived_type;
  typedef std::tuple<std::shared_ptr<LocalfunctionType>> LocalfunctionTupleType;
}; // class LocalGodunovNumericalBoundaryFluxTraits


} // namespace internal


template <class AnalyticalFluxImp, size_t domainDim = AnalyticalFluxImp::dimDomain>
class LocalGodunovNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                   domainDim>>
{
public:
  typedef internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp, domainDim> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange      = Traits::dimRange;

  explicit LocalGodunovNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux, const bool is_linear = false,
                                             const bool reinit_jacobians = true)
    : analytical_flux_(analytical_flux)
    , is_linear_(is_linear)
  {
    if (is_linear_ && (!jacobians_constructed_ || reinit_jacobians))
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return LocalfunctionTupleType();
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& /*local_functions_tuple_entity*/,
                     const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_entity,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_neighbor,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));
    // get flux value
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;
    // get unit outer normal
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif // NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
#ifndef NDEBUG
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
        ++num_zeros;
      else
        DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
#endif // NDEBUG
    }
    assert(num_zeros == dimDomain - 1);
    // calculate return vector
    RangeType ret;
    RangeFieldType vol_intersection = intersection.geometry().volume();
    if (n_ij[coord] > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_[coord].mv(delta_u, negative_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[kk]      = (f_u_i[kk][coord] - negative_waves[kk] * n_ij[coord]) * vol_intersection;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_[coord].mv(delta_u, positive_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[kk]      = (-f_u_i[kk][coord] - positive_waves[kk] * n_ij[coord]) * vol_intersection;
    }
    return ret;
  } // RangeType evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector<EigenMatrixType> jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(
          DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian[ii], 15), dimRange, dimRange));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i, const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector<EigenMatrixType> jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(
          DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian[ii], 15), dimRange, dimRange));
    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector<EigenMatrixType>&& jacobian) const
  {
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver<typename Stuff::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver(
          jacobian[ii].backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues  = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast<size_t>(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      }
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real() * diag_jacobian_neg_tmp.backend()
                                         * eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real() * diag_jacobian_pos_tmp.backend()
                                         * eigenvectors_inverse.real());
      jacobian_neg_[ii] = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_neg_eigen, 15));
      jacobian_pos_[ii] = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_pos_eigen, 15));
    }
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class LocalGodunovNumericalCouplingFlux

template <class AnalyticalFluxImp, size_t dimDomain>
typename internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp, dimDomain>::FluxJacobianRangeType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, dimDomain>::jacobian_neg_(0);

template <class AnalyticalFluxImp, size_t dimDomain>
typename internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp, dimDomain>::FluxJacobianRangeType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, dimDomain>::jacobian_pos_(0);

template <class AnalyticalFluxImp, size_t dimDomain>
bool LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, dimDomain>::jacobians_constructed_(false);

template <class AnalyticalFluxImp>
class LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>
    : public LocalNumericalCouplingFluxInterface<internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp,
                                                                                                   1>>
{
public:
  typedef internal::LocalGodunovNumericalCouplingFluxTraits<AnalyticalFluxImp, 1> Traits;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange  = Traits::dimRange;
  typedef typename DS::Functions::Affine<typename AnalyticalFluxType::FluxDummyEntityType, RangeFieldType, dimRange,
                                         RangeFieldType, dimRange, 1>
      AffineFunctionType;

  explicit LocalGodunovNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux, const bool is_linear,
                                             const bool reinit_jacobians = true)
    : analytical_flux_(analytical_flux)
    , is_linear_(is_linear)
  {
    if (is_linear_ && (!jacobians_constructed_ || reinit_jacobians))
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& /*entity*/) const
  {
    return LocalfunctionTupleType();
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& /*local_functions_tuple_entity*/,
                     const LocalfunctionTupleType& /*local_functions_tuple_neighbor*/,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_entity,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_neighbor,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    // get function values
    const RangeType u_i = local_source_entity.evaluate(intersection.geometryInInside().global(x_intersection));
    const RangeType u_j = local_source_neighbor.evaluate(intersection.geometryInOutside().global(x_intersection));

    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get unit outer normal
    const RangeFieldType n_ij = intersection.unitOuterNormal(x_intersection)[0];
    // calculate return vector
    RangeType ret;
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
    if (n_ij > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_function_.evaluate(u_j - u_i, negative_waves);
      ret = f_u_i + negative_waves;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_function_.evaluate(u_i - u_j, positive_waves);
      ret = positive_waves - f_u_i;
    }
    return ret;
  } // void evaluate(...) const

private:
  void initialize_jacobians() const
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    EigenMatrixType jacobian_eigen(DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian, 15)));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i, const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    EigenMatrixType jacobian_eigen = DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian, 15));

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(EigenMatrixType&& jacobian) const
  {
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    diag_jacobian_pos_tmp.scal(RangeFieldType(0));
    diag_jacobian_neg_tmp.scal(RangeFieldType(0));
    // create EigenSolver
    ::Eigen::EigenSolver<typename Stuff::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver(
        jacobian.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues  = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
    assert(boost::numeric_cast<size_t>(eigenvalues.size()) == dimRange);
    for (size_t jj = 0; jj < dimRange; ++jj) {
      // assert this is real
      assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
      const RangeFieldType eigenvalue = eigenvalues[jj].real();
      if (eigenvalue < 0)
        diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
      else
        diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real() * diag_jacobian_neg_tmp.backend()
                                         * eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real() * diag_jacobian_pos_tmp.backend()
                                         * eigenvectors_inverse.real());
      // set jacobian_neg_ and jacobian_pos_
      jacobian_neg_ = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_neg_eigen, 15));
      jacobian_pos_ = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_pos_eigen, 15));
      jacobian_neg_function_ = AffineFunctionType(jacobian_neg_, RangeType(0), true);
      jacobian_pos_function_ = AffineFunctionType(jacobian_pos_, RangeType(0), true);
    }
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  thread_local static FluxJacobianRangeType jacobian_neg_;
  thread_local static FluxJacobianRangeType jacobian_pos_;
  thread_local static AffineFunctionType jacobian_neg_function_;
  thread_local static AffineFunctionType jacobian_pos_function_;
  thread_local static bool jacobians_constructed_;
  const bool is_linear_;
}; // class LocalGodunovNumericalCouplingFlux< ..., 1 >

template <class AnalyticalFluxImp>
thread_local typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::jacobian_neg_{
        typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType()};

template <class AnalyticalFluxImp>
thread_local typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::jacobian_pos_{
        typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType()};

template <class AnalyticalFluxImp>
thread_local typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::AffineFunctionType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::jacobian_neg_function_(
        LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::AffineFunctionType(
            LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType(0)));

template <class AnalyticalFluxImp>
thread_local typename LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::AffineFunctionType
    LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::jacobian_pos_function_(
        LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::AffineFunctionType(
            LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::FluxJacobianRangeType(0)));

template <class AnalyticalFluxImp>
thread_local bool LocalGodunovNumericalCouplingFlux<AnalyticalFluxImp, 1>::jacobians_constructed_(false);

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp,
          size_t domainDim = AnalyticalBoundaryFluxImp::dimDomain>
class LocalGodunovNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp,
                                                                                             BoundaryValueFunctionImp>>
{
public:
  typedef internal::LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                            domainDim>
      Traits;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange  = Traits::dimRange;

  explicit LocalGodunovNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueFunctionType& boundary_values,
                                             const bool is_linear = false, const bool reinit_jacobians = true)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (is_linear_ && (!jacobians_constructed_ || reinit_jacobians))
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(boundary_values_.local_function(entity));
  }


  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_entity,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
    const RangeType u_i                     = local_source_entity.evaluate(x_intersection_entity_coords);
    const RangeType u_j = std::get<0>(local_functions_tuple)->evaluate(x_intersection_entity_coords);
    // get flux values
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;
    // get unit outer normal
    const auto n_ij = intersection.unitOuterNormal(x_intersection);
    // find direction of unit outer normal
    size_t coord = 0;
#ifndef NDEBUG
    size_t num_zeros = 0;
#endif // NDEBUG
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(1)) || DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(-1)))
        coord = ii;
#ifndef NDEBUG
      else if (DSC::FloatCmp::eq(n_ij[ii], RangeFieldType(0)))
        ++num_zeros;
      else
        DUNE_THROW(Dune::NotImplemented, "Godunov flux is only implemented for axis parallel cube grids");
#endif // NDEBUG
    }
    assert(num_zeros == dimDomain - 1);
    // calculate return vector
    RangeType ret;
    RangeFieldType vol_intersection = intersection.geometry().volume();
    if (n_ij[coord] > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_[coord].mv(delta_u, negative_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[kk]      = (f_u_i[kk][coord] - negative_waves[kk] * n_ij[coord]) * vol_intersection;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_[coord].mv(delta_u, positive_waves);
      for (size_t kk = 0; kk < dimRange; ++kk)
        ret[kk]      = (-f_u_i[kk][coord] - positive_waves[kk] * n_ij[coord]) * vol_intersection;
    }
    return ret;
  } // void evaluate(...) const

private:
  void initialize_jacobians()
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    std::vector<EigenMatrixType> jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian[ii])));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i, const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    std::vector<EigenMatrixType> jacobian_eigen;
    for (size_t ii = 0; ii < dimDomain; ++ii)
      jacobian_eigen.emplace_back(DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian[ii])));
    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(std::vector<EigenMatrixType>&& jacobian) const
  {
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    for (size_t ii = 0; ii < dimDomain; ++ii) {
      diag_jacobian_pos_tmp.scal(RangeFieldType(0));
      diag_jacobian_neg_tmp.scal(RangeFieldType(0));
      // create EigenSolver
      ::Eigen::EigenSolver<typename Stuff::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver(
          jacobian[ii].backend());
      assert(eigen_solver.info() == ::Eigen::Success);
      const auto eigenvalues  = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
      const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
      assert(boost::numeric_cast<size_t>(eigenvalues.size()) == dimRange);
      for (size_t jj = 0; jj < dimRange; ++jj) {
        // assert this is real
        assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
        const RangeFieldType eigenvalue = eigenvalues[jj].real();
        if (eigenvalue < 0)
          diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
        else
          diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      }
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real() * diag_jacobian_neg_tmp.backend()
                                         * eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real() * diag_jacobian_pos_tmp.backend()
                                         * eigenvectors_inverse.real());
      jacobian_neg_[ii] =
          DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(DSC::to_string(jacobian_neg_eigen));
      jacobian_pos_[ii] =
          DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(DSC::to_string(jacobian_pos_eigen));
    }
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  static FluxJacobianRangeType jacobian_neg_;
  static FluxJacobianRangeType jacobian_pos_;
  static bool jacobians_constructed_;
  const bool is_linear_;
}; // class LocalGodunovNumericalBoundaryFlux

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim>
typename internal::LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                           domainDim>::FluxJacobianRangeType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim>::jacobian_neg_(0);

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim>
typename internal::LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                           domainDim>::FluxJacobianRangeType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, domainDim>::jacobian_pos_(0);

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp, size_t domainDim>
bool LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                       domainDim>::jacobians_constructed_(false);

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
class LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp,
                                                                                             BoundaryValueFunctionImp>>
{
public:
  typedef internal::LocalGodunovNumericalBoundaryFluxTraits<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>
      Traits;
  typedef typename Traits::BoundaryValueFunctionType BoundaryValueFunctionType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::FluxRangeType FluxRangeType;
  typedef typename Traits::FluxJacobianRangeType FluxJacobianRangeType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::EigenMatrixType EigenMatrixType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange  = Traits::dimRange;
  typedef typename DS::Functions::Affine<typename AnalyticalFluxType::FluxDummyEntityType, RangeFieldType, dimRange,
                                         RangeFieldType, dimRange, 1>
      AffineFunctionType;
  explicit LocalGodunovNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueFunctionType& boundary_values,
                                             const bool is_linear = false, const bool reinit_jacobians = true)
    : analytical_flux_(analytical_flux)
    , boundary_values_(boundary_values)
    , is_linear_(is_linear)
  {
    if (is_linear_ && (!jacobians_constructed_ || reinit_jacobians))
      initialize_jacobians();
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(boundary_values_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                     const Stuff::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType,
                                                         dimRange, 1>& local_source_entity,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_intersection) const
  {
    const auto x_intersection_entity_coords = intersection.geometryInInside().global(x_intersection);
    const RangeType u_i                     = local_source_entity.evaluate(x_intersection_entity_coords);
    const RangeType u_j = std::get<0>(local_functions_tuple)->evaluate(x_intersection_entity_coords);

    if (!jacobians_constructed_ && is_linear_)
      initialize_jacobians();

    if (!is_linear_) { // use simple linearized Riemann solver, LeVeque p.316
      reinitialize_jacobians(u_i, u_j);
    }
    // get unit outer normal
    const RangeFieldType n_ij = intersection.unitOuterNormal(x_intersection)[0];
    // calculate return vector
    RangeType ret;
    const FluxRangeType f_u_i = analytical_flux_.evaluate(u_i);
    if (n_ij > 0) {
      RangeType negative_waves(RangeFieldType(0));
      jacobian_neg_function_.evaluate(u_j - u_i, negative_waves);
      ret = f_u_i + negative_waves;
    } else {
      RangeType positive_waves(RangeFieldType(0));
      jacobian_pos_function_.evaluate(u_i - u_j, positive_waves);
      ret = positive_waves - f_u_i;
    }
    return ret;
  } // void evaluate(...) const

private:
  void initialize_jacobians() const
  {
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(RangeType(0)));
    EigenMatrixType jacobian_eigen(DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian, 15)));
    calculate_jacobians(std::move(jacobian_eigen));
    jacobians_constructed_ = true;
  } // void initialize_jacobians()

  void reinitialize_jacobians(const RangeType& u_i, const RangeType& u_j) const
  {
    // calculate jacobian as jacobian(0.5*(u_i+u_j)
    RangeType u_mean = u_i + u_j;
    u_mean *= RangeFieldType(0.5);
    const FluxJacobianRangeType jacobian(analytical_flux_.jacobian(u_mean));
    EigenMatrixType jacobian_eigen = DSC::from_string<EigenMatrixType>(DSC::to_string(jacobian, 15));

    // calculate jacobian_neg and jacobian_pos
    calculate_jacobians(std::move(jacobian_eigen));
  } // void reinitialize_jacobians(...)

  void calculate_jacobians(EigenMatrixType&& jacobian) const
  {
    EigenMatrixType diag_jacobian_pos_tmp(dimRange, dimRange, RangeFieldType(0));
    EigenMatrixType diag_jacobian_neg_tmp(dimRange, dimRange, RangeFieldType(0));
    diag_jacobian_pos_tmp.scal(RangeFieldType(0));
    diag_jacobian_neg_tmp.scal(RangeFieldType(0));
    // create EigenSolver
    ::Eigen::EigenSolver<typename Stuff::LA::EigenDenseMatrix<RangeFieldType>::BackendType> eigen_solver(
        jacobian.backend());
    assert(eigen_solver.info() == ::Eigen::Success);
    const auto eigenvalues  = eigen_solver.eigenvalues(); // <- this should be an Eigen vector of std::complex
    const auto eigenvectors = eigen_solver.eigenvectors(); // <- this should be an Eigen matrix of std::complex
    assert(boost::numeric_cast<size_t>(eigenvalues.size()) == dimRange);
    for (size_t jj = 0; jj < dimRange; ++jj) {
      // assert this is real
      assert(std::abs(eigenvalues[jj].imag()) < 1e-15);
      const RangeFieldType eigenvalue = eigenvalues[jj].real();
      if (eigenvalue < 0)
        diag_jacobian_neg_tmp.set_entry(jj, jj, eigenvalue);
      else
        diag_jacobian_pos_tmp.set_entry(jj, jj, eigenvalue);
      const auto eigenvectors_inverse = eigenvectors.inverse();
      EigenMatrixType jacobian_neg_eigen(eigenvectors.real() * diag_jacobian_neg_tmp.backend()
                                         * eigenvectors_inverse.real());
      EigenMatrixType jacobian_pos_eigen(eigenvectors.real() * diag_jacobian_pos_tmp.backend()
                                         * eigenvectors_inverse.real());
      // set jacobian_neg_ and jacobian_pos_
      jacobian_neg_ = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_neg_eigen, 15));
      jacobian_pos_ = DSC::from_string<Dune::FieldMatrix<RangeFieldType, dimRange, dimRange>>(
          DSC::to_string(jacobian_pos_eigen, 15));
      jacobian_neg_function_ = AffineFunctionType(jacobian_neg_, RangeType(0), true);
      jacobian_pos_function_ = AffineFunctionType(jacobian_pos_, RangeType(0), true);
    }
  } // void calculate_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  const BoundaryValueFunctionType& boundary_values_;
  thread_local static FluxJacobianRangeType jacobian_neg_;
  thread_local static FluxJacobianRangeType jacobian_pos_;
  thread_local static AffineFunctionType jacobian_neg_function_;
  thread_local static AffineFunctionType jacobian_pos_function_;
  thread_local static bool jacobians_constructed_;
  const bool is_linear_;
}; // class LocalGodunovNumericalBoundaryFlux< ..., 1 >

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
thread_local typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                        1>::FluxJacobianRangeType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::jacobian_neg_{
        typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                   1>::FluxJacobianRangeType()};

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
thread_local typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                        1>::FluxJacobianRangeType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::jacobian_pos_{
        typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                   1>::FluxJacobianRangeType()};

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
thread_local typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                        1>::AffineFunctionType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::jacobian_neg_function_(
        LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::AffineFunctionType(
            LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                              1>::FluxJacobianRangeType(0)));

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
thread_local typename LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                                        1>::AffineFunctionType
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::jacobian_pos_function_(
        LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::AffineFunctionType(
            LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp,
                                              1>::FluxJacobianRangeType(0)));

template <class AnalyticalBoundaryFluxImp, class BoundaryValueFunctionImp>
thread_local bool
    LocalGodunovNumericalBoundaryFlux<AnalyticalBoundaryFluxImp, BoundaryValueFunctionImp, 1>::jacobians_constructed_(
        false);


#else // HAVE_EIGEN

template <class AnalyticalFluxImp, size_t domainDim>
class LocalGodunovNumericalCouplingFlux
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

template <class AnalyticalFluxImp, class BoundaryValueFunctionType, size_t domainDim>
class LocalGodunovNumericalBoundaryFlux
{
  static_assert(AlwaysFalse<AnalyticalFluxImp>::value, "You are missing eigen!");
};

#endif // HAVE_EIGEN


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
