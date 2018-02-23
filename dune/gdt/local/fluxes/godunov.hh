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

#include "interfaces.hh"
#include "base.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp, class Traits>
class GodunovLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueType, class BoundaryInfoType, class Traits>
class GodunovLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp>
class GodunovLocalNumericalCouplingFluxTraits : public NumericalCouplingFluxTraitsBase<AnalyticalFluxImp>
{
  typedef NumericalCouplingFluxTraitsBase<AnalyticalFluxImp> BaseType;

public:
  using typename BaseType::AnalyticalFluxLocalfunctionType;
  typedef GodunovLocalNumericalCouplingFlux<AnalyticalFluxImp, GodunovLocalNumericalCouplingFluxTraits> derived_type;
  typedef std::tuple<std::shared_ptr<AnalyticalFluxLocalfunctionType>> LocalfunctionTupleType;
}; // class GodunovLocalNumericalCouplingFluxTraits

template <class AnalyticalFluxImp, class BoundaryValueImp, class BoundaryInfoImp>
class GodunovLocalDirichletNumericalBoundaryFluxTraits
    : public NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp, BoundaryInfoImp>
{
  typedef NumericalBoundaryFluxTraitsBase<AnalyticalFluxImp, BoundaryValueImp, BoundaryInfoImp> BaseType;

public:
  typedef GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxImp,
                                                     BoundaryValueImp,
                                                     BoundaryInfoImp,
                                                     GodunovLocalDirichletNumericalBoundaryFluxTraits>
      derived_type;
  using typename BaseType::AnalyticalFluxLocalfunctionType;
  using typename BaseType::BoundaryValueLocalfunctionType;
  typedef std::tuple<std::shared_ptr<AnalyticalFluxLocalfunctionType>, std::shared_ptr<BoundaryValueLocalfunctionType>>
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
  typedef typename Traits::AnalyticalFluxLocalfunctionType AnalyticalFluxLocalfunctionType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;
  typedef typename Dune::FieldMatrix<RangeFieldType, dimRange, dimRange> MatrixType;
  typedef typename XT::LA::CommonSparseMatrix<RangeFieldType> SparseMatrixType;
  typedef typename XT::LA::CommonSparseMatrixCsc<RangeFieldType> CscSparseMatrixType;
  typedef typename XT::LA::EigenSolver<MatrixType> EigenSolverType;
  typedef typename XT::LA::EigenSolverOptions<MatrixType> EigenSolverOptionsType;
  typedef typename XT::LA::MatrixInverterOptions<MatrixType> MatrixInverterOptionsType;
  typedef typename AnalyticalFluxLocalfunctionType::StateRangeType StateRangeType;
  typedef typename Dune::FieldVector<MatrixType, dimDomain> JacobianRangeType;

  typedef FieldVector<CscSparseMatrixType, dimDomain> JacobiansType;

  explicit GodunovFluxImplementation(const AnalyticalFluxType& analytical_flux,
                                     XT::Common::Parameter param,
                                     const bool is_linear = false,
                                     const bool boundary = false)
    : analytical_flux_(analytical_flux)
    , param_inside_(param)
    , param_outside_(param)
    , is_linear_(is_linear)
  {
    param_inside_.set("boundary", {0.}, true);
    param_outside_.set("boundary", {double(boundary)}, true);
    if (is_instantiated_)
      DUNE_THROW(InvalidStateException,
                 "This class uses several static variables to save its state between time "
                 "steps, so using several instances at the same time may result in undefined "
                 "behavior!");
    is_instantiated_ = true;
  }

  ~GodunovFluxImplementation()
  {
    is_instantiated_ = false;
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
    initialize_jacobians(direction, local_functions_tuple, x_in_inside_coords, u_i, u_j);

    // get jump at the intersection
    const RangeType delta_u = u_i - u_j;
    // calculate waves
    RangeType waves(0);
    n_ij[direction] > 0 ? jacobian_neg()[direction].mtv(delta_u, waves) : jacobian_pos()[direction].mtv(delta_u, waves);
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

  static void reset()
  {
    ++initialization_count_;
  }

private:
  // solve R^T F = P^T D A^T where is P is the permutation matrix and D is the diagonal eigenvalue matrix
  static void solve_upper_triangular_transposed(CscSparseMatrixType& F,
                                                const SparseMatrixType& R,
                                                const FieldVector<int, dimRange>& permutations,
                                                const StateRangeType& eigenvalues,
                                                const SparseMatrixType& A)
  {
    F.clear();
    StateRangeType rhs(0.), x(0.);
    for (size_t cc = 0; cc < dimRange; ++cc) {
      // apply permutations and scaling to column of A^T and store in vector
      std::fill(rhs.begin(), rhs.end(), 0.);
      for (size_t kk = A.outer_index_ptr()[cc]; kk < A.outer_index_ptr()[cc + 1]; ++kk) {
        size_t index = A.inner_index_ptr()[kk];
        auto eigval = eigenvalues[index];
        if (XT::Common::FloatCmp::ne(eigval, 0.))
          rhs[permutations[index]] = A.entries()[kk] * eigval;
      } // kk
      XT::LA::solve_upper_triangular(R, x, rhs);
      F.start_column();
      for (size_t rr = 0; rr < dimRange; ++rr)
        if (XT::Common::FloatCmp::ne(x[rr], 0.))
          F.push_entry(rr, x[rr]);
      F.end_column();
    } // cc
  } // void solve_upper_triangular_transposed(...)

  static XT::Common::Configuration create_eigensolver_options()
  {
    XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options(EigenSolverOptionsType::types()[0]);
    // XT::Common::Configuration eigensolver_options = EigenSolverOptionsType::options("shifted_qr");
    eigensolver_options["assert_eigendecomposition"] = "1e-6";
    eigensolver_options["assert_real_eigendecomposition"] = "1e-6";
    eigensolver_options["disable_checks"] = "true";
    XT::Common::Configuration matrix_inverter_options = MatrixInverterOptionsType::options();
    matrix_inverter_options["post_check_is_left_inverse"] = "1e-6";
    matrix_inverter_options["post_check_is_right_inverse"] = "1e-6";
    eigensolver_options.add(matrix_inverter_options, "matrix-inverter");
    return eigensolver_options;
  } // ... create_eigensolver_options()

  // use simple linearized Riemann solver, LeVeque p.316
  void initialize_jacobians(const size_t direction,
                            const LocalfunctionTupleType& local_functions_tuple,
                            const DomainType& x_local,
                            const RangeType& u_i,
                            const RangeType& u_j) const
  {
    if (local_initialization_counts_[direction] != initialization_count_) {
      // calculate jacobian as jacobian(0.5*(u_i+u_j)
      RangeType u_mean = u_i + u_j;
      u_mean *= RangeFieldType(0.5);
      const auto& local_flux = std::get<0>(local_functions_tuple);
      thread_local std::unique_ptr<MatrixType> jacobian;
      thread_local RangeType tau;
      if (!jacobian) {
        jacobian = XT::Common::make_unique<MatrixType>();
      }
      helper<dimDomain>::get_jacobian(direction, local_flux, x_local, u_mean, *jacobian, param_inside_);
      // get matrix of eigenvectors A and eigenvalues
      static auto eigensolver_options = create_eigensolver_options();
      const auto eigen_solver = EigenSolverType(*jacobian, eigensolver_options);
      const auto eigenvalues = eigen_solver.real_eigenvalues();
      StateRangeType eigvals_neg(0.);
      StateRangeType eigvals_pos(0.);
      for (size_t ii = 0; ii < eigenvalues.size(); ++ii) {
        if (eigenvalues[ii] < 0.)
          eigvals_neg[ii] = eigenvalues[ii];
        else
          eigvals_pos[ii] = eigenvalues[ii];
      }
      thread_local std::unique_ptr<MatrixType> eigenvectors_dense = XT::Common::make_unique<MatrixType>();
      *eigenvectors_dense = eigen_solver.real_eigenvectors();
      thread_local SparseMatrixType eigenvectors(dimRange, dimRange, size_t(0));
      eigenvectors = *eigenvectors_dense;
      // calculate QR decomposition with column pivoting A = QRP^T
      FieldVector<int, dimRange> permutations;
      XT::LA::qr(*eigenvectors_dense, tau, permutations);
      static auto upper_triangular_pattern =
          XT::Common::triangular_pattern(dimRange, dimRange, XT::Common::MatrixPattern::upper_triangular);
      thread_local SparseMatrixType R(dimRange, dimRange, size_t(0));
      R.assign(*eigenvectors_dense, upper_triangular_pattern);
      XT::LA::calculate_q_from_qr(*eigenvectors_dense, tau);
      const auto& Q = *eigenvectors_dense;
      FieldVector<int, dimRange> inverse_permutations;
      for (size_t ii = 0; ii < dimRange; ++ii)
        inverse_permutations[permutations[ii]] = ii;
      // we want to compute C_{+,-} = A D_{+,-} A^{-1}, where D_+ is the diagonal matrix containing only positive
      // eigenvalues and D_- contains only negative eigenvalues. Using the QR decomposition, we get
      // C = A D A^{-1} = A D (QRP^T)^{-1} = A D P R^{-1} Q^T. For the transpose of C, we get
      // C^T =  Q R^{-T} P^T D A^T. We calculate C^T by first solving (column-wise)
      // R^T F = P^T D A^T and then computing C^T = Q F;
      thread_local CscSparseMatrixType F_neg(dimRange, dimRange, size_t(0));
      solve_upper_triangular_transposed(F_neg, R, inverse_permutations, eigvals_neg, eigenvectors);
      jacobian_neg()[direction] = Q;
      jacobian_pos()[direction].deep_copy(jacobian_neg()[direction]);
      jacobian_neg()[direction].rightmultiply(F_neg);
      auto& F_pos = F_neg;
      solve_upper_triangular_transposed(F_pos, R, inverse_permutations, eigvals_pos, eigenvectors);
      jacobian_pos()[direction].rightmultiply(F_pos);
      if (is_linear_) {
        jacobian = nullptr;
        ++local_initialization_counts_[direction];
      }
    } // if (local_initialization_counts_[direction] != initialization_count_)
  } // void calculate_jacobians(...)

  template <size_t domainDim = dimDomain, class anything = void>
  struct helper
  {
    static void get_jacobian(const size_t direction,
                             const std::shared_ptr<AnalyticalFluxLocalfunctionType>& local_func,
                             const DomainType& x_in_inside_coords,
                             const StateRangeType& u,
                             MatrixType& ret,
                             const XT::Common::Parameter& param)
    {
      local_func->partial_u_col(direction, x_in_inside_coords, u, ret, param);
    }
  };

  template <class anything>
  struct helper<1, anything>
  {
    static void get_jacobian(const size_t direction,
                             const std::shared_ptr<AnalyticalFluxLocalfunctionType>& local_func,
                             const DomainType& x_in_inside_coords,
                             const StateRangeType& u,
                             MatrixType& ret,
                             const XT::Common::Parameter& param)
    {
      assert(direction == 0);
      local_func->partial_u(x_in_inside_coords, u, ret, param);
    }
  };

  const AnalyticalFluxType& analytical_flux_;
  XT::Common::Parameter param_inside_;
  XT::Common::Parameter param_outside_;

  // work around gcc bug 66944
  static JacobiansType& jacobian_neg()
  {
    static thread_local JacobiansType jacobian_neg_(CscSparseMatrixType(dimRange, dimRange, size_t(0)));
    return jacobian_neg_;
  }

  static JacobiansType& jacobian_pos()
  {
    static thread_local JacobiansType jacobian_pos_(CscSparseMatrixType(dimRange, dimRange, size_t(0)));
    return jacobian_pos_;
  }

  //  static thread_local JacobiansType jacobian_neg_;
  //  static thread_local JacobiansType jacobian_pos_;

  const bool is_linear_;
  static bool is_instantiated_;
  static std::atomic<size_t> initialization_count_;
  static thread_local FieldVector<size_t, dimDomain> local_initialization_counts_;
}; // class GodunovFluxImplementation

// template <class Traits>
// thread_local
//    typename GodunovFluxImplementation<Traits>::JacobiansType GodunovFluxImplementation<Traits>::jacobian_neg_;

// template <class Traits>
// thread_local
//    typename GodunovFluxImplementation<Traits>::JacobiansType GodunovFluxImplementation<Traits>::jacobian_pos_;

template <class Traits>
bool GodunovFluxImplementation<Traits>::is_instantiated_(false);

template <class Traits>
std::atomic<size_t> GodunovFluxImplementation<Traits>::initialization_count_(1);

template <class Traits>
thread_local FieldVector<size_t, GodunovFluxImplementation<Traits>::dimDomain>
    GodunovFluxImplementation<Traits>::local_initialization_counts_(0);


} // namespace internal


template <class AnalyticalFluxImp, class Traits = internal::GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp>>
class GodunovLocalNumericalCouplingFlux : public LocalNumericalCouplingFluxInterface<Traits>
{
public:
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  static constexpr size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit GodunovLocalNumericalCouplingFlux(const AnalyticalFluxType& analytical_flux,
                                             const XT::Common::Parameter& param,
                                             const bool is_linear = false)
    : implementation_(analytical_flux, param, is_linear, false)
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
          class BoundaryInfoImp,
          class Traits = internal::
              GodunovLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp, BoundaryValueImp, BoundaryInfoImp>>
class GodunovLocalDirichletNumericalBoundaryFlux : public LocalNumericalBoundaryFluxInterface<Traits>
{
  typedef LocalNumericalBoundaryFluxInterface<Traits> InterfaceType;

public:
  typedef typename Traits::BoundaryValueType BoundaryValueType;
  typedef typename Traits::BoundaryInfoType BoundaryInfoType;
  typedef typename Traits::LocalfunctionTupleType LocalfunctionTupleType;
  typedef typename Traits::EntityType EntityType;
  typedef typename Traits::DomainFieldType DomainFieldType;
  typedef typename Traits::RangeFieldType RangeFieldType;
  typedef typename Traits::AnalyticalFluxType AnalyticalFluxType;
  typedef typename Traits::RangeType RangeType;
  typedef typename Traits::DomainType DomainType;
  static const size_t dimDomain = Traits::dimDomain;
  static const size_t dimRange = Traits::dimRange;

  explicit GodunovLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                      const BoundaryValueType& boundary_values,
                                                      const BoundaryInfoType& boundary_info,
                                                      const XT::Common::Parameter& param,
                                                      const bool is_linear = false)
    : InterfaceType(boundary_values, boundary_info)
    , implementation_(analytical_flux, param, is_linear, true)
  {
  }

  LocalfunctionTupleType local_functions(const EntityType& entity) const
  {
    return std::make_tuple(implementation_.analytical_flux().local_function(entity),
                           std::get<0>(boundary_values_)->local_function(entity));
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
