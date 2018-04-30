// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2016 - 2017)
//   Rene Milk       (2016 - 2018)
//   Tobias Leibner  (2016 - 2017)

#ifndef DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
#define DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH

#include <tuple>
#include <memory>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include <dune/gdt/operators/fv/reconstruction/linear.hh>

#include "interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class Traits>
class GodunovLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueType, class Traits>
class GodunovLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class MatrixType, class VectorType, size_t dimRange, size_t num_jacobians>
class GodunovJacobianWrapper : public JacobianWrapper<MatrixType, VectorType, dimRange, num_jacobians>
{
  using BaseType = JacobianWrapper<MatrixType, VectorType, dimRange, num_jacobians>;
  using typename BaseType::V;

public:
  GodunovJacobianWrapper()
    : eigvals_neg_(V::create(dimRange))
    , eigvals_pos_(V::create(dimRange))
  {
  }

  void compute()
  {
    for (size_t dd = 0; dd < num_jacobians; ++dd)
      compute(dd);
  }

  void compute(const size_t dd)
  {
    BaseType::compute(dd);
    std::fill(eigvals_neg_[dd].begin(), eigvals_neg_[dd].end(), 0.);
    std::fill(eigvals_pos_[dd].begin(), eigvals_pos_[dd].end(), 0.);
    for (size_t ii = 0; ii < eigenvalues_[dd].size(); ++ii)
      (eigenvalues_[dd][ii] < 0 ? eigvals_neg_[dd][ii] : eigvals_pos_[dd][ii]) = eigenvalues_[dd][ii];
  }

  const VectorType& negative_eigenvalues(const size_t direction) const
  {
    return eigvals_neg_[direction];
  }

  const VectorType& positive_eigenvalues(const size_t direction) const
  {
    return eigvals_pos_[direction];
  }

protected:
  using BaseType::eigenvalues_;
  using BaseType::computed_;
  FieldVector<VectorType, num_jacobians> eigvals_neg_;
  FieldVector<VectorType, num_jacobians> eigvals_pos_;
};


template <class AnalyticalFluxImp>
class GodunovLocalNumericalCouplingFluxTraits : public NumericalCouplingFluxTraitsBase<AnalyticalFluxImp>
{
  typedef NumericalCouplingFluxTraitsBase<AnalyticalFluxImp> BaseType;

public:
  using typename BaseType::AnalyticalFluxLocalfunctionType;
  typedef GodunovLocalNumericalCouplingFlux<AnalyticalFluxImp, GodunovLocalNumericalCouplingFluxTraits> derived_type;
  typedef std::tuple<std::unique_ptr<AnalyticalFluxLocalfunctionType>> LocalfunctionTupleType;
}; // class GodunovLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp>
class GodunovLocalDirichletNumericalBoundaryFluxTraits
    : public NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp>
{
  typedef NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp> BaseType;

public:
  typedef GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                     BoundaryValueImp,
                                                     GodunovLocalDirichletNumericalBoundaryFluxTraits>
      derived_type;
  using typename BaseType::AnalyticalFluxLocalfunctionType;
  using typename BaseType::LocalBoundaryValueType;
  typedef std::tuple<std::unique_ptr<AnalyticalFluxLocalfunctionType>, std::unique_ptr<LocalBoundaryValueType>>
      LocalfunctionTupleType;
}; // class GodunovLocalDirichletNumericalBoundaryFluxTraits

template <class Traits>
class GodunovFluxImplementation
{
public:
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  using MatrixType = FieldMatrix<DomainFieldType, dimRange, dimRange>;
  using VectorType = std::vector<RangeFieldType>;
  using JacobianWrapperType = GodunovJacobianWrapper<MatrixType, VectorType, dimRange, dimDomain>;

  explicit GodunovFluxImplementation(const AnalyticalFluxType& analytical_flux,
                                     XT::Common::Parameter param,
                                     XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper,
                                     const bool boundary = false)
    : analytical_flux_(analytical_flux)
    , param_inside_(param)
    , param_outside_(param)
    , jacobian_wrapper_(jacobian_wrapper)
  {
    param_inside_.set("boundary", {0.}, true);
    param_outside_.set("boundary", {double(boundary)}, true);
  }

  template <class IntersectionType>
  RangeType evaluate(const LocalfunctionTupleType& local_functions_tuple,
                     const IntersectionType& intersection,
                     const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords,
                     const DomainType& x_in_inside_coords,
                     const DomainType& /*x_in_outside_coords*/,
                     const RangeType& u_i,
                     const RangeType& u_j) const
  {
    size_t direction = intersection.indexInInside() / 2;
    // get unit outer normal
    const auto n_ij = intersection.unitOuterNormal(x_in_intersection_coords);
    assert(XT::Common::FloatCmp::eq(std::abs(n_ij[direction]), 1.));

    // intialize jacobians
    auto& jac = *jacobian_wrapper_;
    initialize_jacobians(jac, direction, local_functions_tuple, x_in_inside_coords, u_i, u_j);

    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;

    // calculate waves

    // apply A D_{+,-} A^{-1}, where A is the matrix of eigenvectors and D_{+/-} is the diagonal matrix containing only
    // positive/negative eigenvalues of A
    RangeType waves, tmp_vec;
    jac.apply_inverse_eigenvectors(direction, delta_u, tmp_vec);
    const auto& eigenvalues =
        n_ij[direction] > 0 ? jac.negative_eigenvalues(direction) : jac.positive_eigenvalues(direction);
    for (size_t ii = 0; ii < dimRange; ++ii)
      tmp_vec[ii] *= eigenvalues[ii];
    jac.apply_eigenvectors(direction, tmp_vec, waves);

    // calculate flux
    const auto& local_flux = std::get<0>(local_functions_tuple);
    RangeType ret = local_flux->evaluate_col(direction, x_in_inside_coords, u_i, param_inside_);
    ret -= waves;
    ret *= n_ij[direction];
    return ret;
  } // RangeType evaluate(...) const

  const AnalyticalFluxType& analytical_flux() const
  {
    return analytical_flux_;
  }

private:
  // use simple linearized Riemann solver, LeVeque p.316
  void initialize_jacobians(JacobianWrapperType& jac,
                            const size_t direction,
                            const LocalfunctionTupleType& local_functions_tuple,
                            const DomainType& x_local,
                            const RangeType& u_i,
                            const RangeType& u_j) const
  {
    // get jacobian
    if (!jac.computed(direction) || !analytical_flux_.is_affine()) {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      auto u_mean = (u_i + u_j) * 0.5;
      const auto& local_flux = std::get<0>(local_functions_tuple);
      local_flux->partial_u_col(direction, x_local, u_mean, jac.jacobian(direction), param_inside_);
      jac.compute(direction);
    }
  } // void initialize_jacobians(...)

  const AnalyticalFluxType& analytical_flux_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;
  XT::Common::PerThreadValue<JacobianWrapperType>& jacobian_wrapper_;
}; // class GodunovFluxImplementation


} // namespace internal


template <class AnalyticalFluxImp,
          class TraitsImp = internal::GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp>>
class GodunovLocalNumericalCouplingFlux : public LocalNumericalCouplingFluxInterface<TraitsImp>
{
public:
  using Traits = TraitsImp;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  template <class JacobianWrapperType>
  GodunovLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                    const XT::Common::Parameter& param,
                                    JacobianWrapperType&& jacobian_wrapper)
    : implementation_(analytical_flux, param, std::forward<JacobianWrapperType>(jacobian_wrapper), false)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(implementation_.analytical_flux().local_function(entity));
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
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const auto x_in_outside_coords = intersection.geometryInOutside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = local_source_neighbor.evaluate(x_in_outside_coords);
    return implementation_.evaluate(local_functions_tuple_entity,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_outside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

  // clear static variables
  static void reset()
  {
    internal::GodunovFluxImplementation<Traits>::reset();
  }

private:
  const internal::GodunovFluxImplementation<Traits> implementation_;
}; // class GodunovLocalNumericalCouplingFlux

/**
*  \brief  Godunov flux evaluation for Dirichlet boundary intersections.
*  \see    GodunovLocalNumericalCouplingFlux
*/
template <class AnalyticalFluxImp,
          class BoundaryValueImp,
          class Traits =
              internal::GodunovLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp, BoundaryValueImp>>
class GodunovLocalDirichletNumericalBoundaryFlux : public LocalNumericalBoundaryFluxInterface<Traits>
{
  typedef LocalNumericalBoundaryFluxInterface<Traits> InterfaceType;

public:
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  template <class JacobianWrapperType>
  GodunovLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                             const BoundaryValueType& boundary_values,
                                             const XT::Common::Parameter& param,
                                             JacobianWrapperType&& jacobian_wrapper)
    : InterfaceType(boundary_values)
    , implementation_(analytical_flux, param, std::forward<JacobianWrapperType>(jacobian_wrapper), true)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(implementation_.analytical_flux().local_function(entity),
                           boundary_values_.local_function(entity));
  }

  template <class IntersectionType>
  RangeType evaluate(
      const LocalfunctionTupleType& local_functions_tuple,
      const XT::Functions::LocalfunctionInterface<EntityType, DomainFieldType, dimDomain, RangeFieldType, dimRange, 1>&
          local_source_entity,
      const IntersectionType& intersection,
      const Dune::FieldVector<DomainFieldType, dimDomain - 1>& x_in_intersection_coords) const

  {
    const auto values = InterfaceType::template get_values<IntersectionType, 1>(
        local_functions_tuple, local_source_entity, intersection, x_in_intersection_coords);
    return implementation_.evaluate(local_functions_tuple,
                                    intersection,
                                    x_in_intersection_coords,
                                    std::get<2>(values),
                                    std::get<2>(values),
                                    std::get<0>(values),
                                    std::get<1>(values));
  } // RangeType evaluate(...) const

  // clear static variables
  static void reset()
  {
    internal::GodunovFluxImplementation<Traits>::reset();
  }

private:
  using InterfaceType::boundary_values_;
  const internal::GodunovFluxImplementation<Traits> implementation_;
}; // class GodunovLocalDirichletNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
