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

#ifndef DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
#define DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH

#include <tuple>
#include <memory>

#include <dune/xt/common/memory.hh>

#include <dune/xt/functions/interfaces.hh>

#include <dune/xt/la/algorithms/qr.hh>
#include <dune/xt/la/eigen-solver.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forwards
template <class AnalyticalFluxImp>
class GodunovLocalNumericalCouplingFlux;

template <class AnalyticalFluxImp, class BoundaryValueType>
class GodunovLocalDirichletNumericalBoundaryFlux;


namespace internal {


template <class AnalyticalFluxImp>
class GodunovLocalNumericalCouplingFluxTraits
{
  //  static_assert(is_analytical_flux<AnalyticalFluxImp>::value,
  //                "AnalyticalFluxImp has to be derived from AnalyticalFluxInterface");

public:
  typedef AnalyticalFluxImp AnalyticalFluxType;
  typedef GodunovLocalNumericalCouplingFlux<AnalyticalFluxType> derived_type;
  typedef typename AnalyticalFluxType::EntityType EntityType;
  typedef typename AnalyticalFluxType::DomainFieldType DomainFieldType;
  typedef typename AnalyticalFluxType::DomainType DomainType;
  typedef typename AnalyticalFluxType::RangeFieldType RangeFieldType;
  typedef typename AnalyticalFluxType::StateRangeType RangeType;
  typedef typename AnalyticalFluxType::LocalfunctionType AnalyticalFluxLocalfunctionType;
  typedef std::tuple<std::shared_ptr<AnalyticalFluxLocalfunctionType>> LocalfunctionTupleType;
  static const size_t dimDomain = AnalyticalFluxType::dimDomain;
  static const size_t dimRange = AnalyticalFluxType::dimRange;
}; // class GodunovLocalNumericalCouplingFluxTraits

template <class AnalyticalBoundaryFluxImp, class BoundaryValueImp>
class GodunovLocalDirichletNumericalBoundaryFluxTraits
    : public GodunovLocalNumericalCouplingFluxTraits<AnalyticalBoundaryFluxImp>
{
  typedef GodunovLocalNumericalCouplingFluxTraits<AnalyticalBoundaryFluxImp> BaseType;

public:
  typedef AnalyticalBoundaryFluxImp AnalyticalFluxType;
  typedef BoundaryValueImp BoundaryValueType;
  typedef typename BoundaryValueType::LocalfunctionType LocalfunctionType;
  typedef GodunovLocalDirichletNumericalBoundaryFlux<AnalyticalFluxType, BoundaryValueType> derived_type;
  using typename BaseType::AnalyticalFluxLocalfunctionType;
  typedef std::tuple<std::shared_ptr<AnalyticalFluxLocalfunctionType>, std::shared_ptr<LocalfunctionType>>
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
                                                const FieldVector<size_t, dimRange>& permutations,
                                                const StateRangeType& eigenvalues,
                                                const SparseMatrixType& A)
  {
    F.clear();
    StateRangeType rhs(0.);
    for (size_t cc = 0; cc < dimRange; ++cc) {
      // apply permutations and scaling to column of A^T and store in vector
      std::fill(rhs.begin(), rhs.end(), 0.);
      size_t end = A.row_pointers()[cc + 1];
      for (size_t kk = A.row_pointers()[cc]; kk < end; ++kk) {
        size_t index = A.column_indices()[kk];
        auto eigval = eigenvalues[index];
        if (XT::Common::FloatCmp::ne(eigval, 0.))
          rhs[permutations[index]] = A.entries()[kk] * eigval;
      } // kk
      // forward solve
      for (size_t ii = 0; ii < dimRange; ++ii) {
        if (XT::Common::FloatCmp::ne(rhs[ii], 0.)) {
          size_t kk = R.row_pointers()[ii];
          const auto& diag_ii = R.entries()[kk];
          const auto x_ii = rhs[ii] / diag_ii;
          F.entries().push_back(x_ii);
          F.row_indices().push_back(ii);
          ++kk;
          for (; kk < R.row_pointers()[ii + 1]; ++kk)
            rhs[R.column_indices()[kk]] -= R.entries()[kk] * x_ii;
        }
      } // ii
      F.column_pointers()[cc + 1] = F.row_indices().size();
    } // cc
  } // void solve_upper_triangular_transposed(...)

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
      thread_local std::unique_ptr<MatrixType> Qdense;
      if (!jacobian) {
        jacobian = XT::Common::make_unique<MatrixType>();
        Qdense = XT::Common::make_unique<MatrixType>();
      }
      helper<dimDomain>::get_jacobian(direction, local_flux, x_local, u_mean, *jacobian, param_inside_);
      // get matrix of eigenvectors A and eigenvalues
      static XT::Common::Configuration eigensolver_options(
          {"type", "check_for_inf_nan", "check_evs_are_real", "check_evs_are_positive", "check_eigenvectors_are_real"},
          {EigenSolverType::types()[0].c_str(), "1", "1", "0", "1"});
      const auto eigen_solver = EigenSolverType(*jacobian);
      const auto eigenvalues = eigen_solver.eigenvalues(eigensolver_options);
      StateRangeType eigvals_neg(0.);
      StateRangeType eigvals_pos(0.);
      for (size_t ii = 0; ii < eigenvalues.size(); ++ii) {
        if (eigenvalues[ii].real() < 0.)
          eigvals_neg[ii] = eigenvalues[ii].real();
        else
          eigvals_pos[ii] = eigenvalues[ii].real();
      }
      const auto eigenvectors_dense = eigen_solver.real_eigenvectors_as_matrix(eigensolver_options);
      thread_local SparseMatrixType eigenvectors(dimRange, dimRange, size_t(0));
      eigenvectors = *eigenvectors_dense;
      // calculate QR decomposition with column pivoting A = QRP^T
      FieldVector<size_t, dimRange> permutations;
      XT::LA::qr_decomposition(*eigenvectors_dense, permutations, *Qdense);
      FieldVector<size_t, dimRange> inverse_permutations;
      for (size_t ii = 0; ii < dimRange; ++ii)
        inverse_permutations[permutations[ii]] = ii;
      thread_local SparseMatrixType R(dimRange, dimRange, size_t(0));
      R = *eigenvectors_dense;
      // we want to compute C_{+,-} = A D_{+,-} A^{-1}, where D_+ is the diagonal matrix containing only positive
      // eigenvalues and D_- contains only negative eigenvalues. Using the QR decomposition, we get
      // C = A D A^{-1} = A D (QRP^T)^{-1} = A D P R^{-1} Q^T. For the transpose of C, we get
      // C^T =  Q R^{-T} P^T D A^T. We calculate C^T by first solving (column-wise)
      // R^T F = P^T D A^T and then computing C^T = Q F;
      thread_local CscSparseMatrixType F_neg(dimRange, dimRange, size_t(0));
      thread_local CscSparseMatrixType F_pos(dimRange, dimRange, size_t(0));
      solve_upper_triangular_transposed(F_neg, R, inverse_permutations, eigvals_neg, eigenvectors);
      solve_upper_triangular_transposed(F_pos, R, inverse_permutations, eigvals_pos, eigenvectors);

      jacobian_neg()[direction] = *Qdense;
      jacobian_pos()[direction].deep_copy(jacobian_neg()[direction]);
      jacobian_neg()[direction].rightmultiply(F_neg);
      jacobian_pos()[direction].rightmultiply(F_pos);

      if (is_linear_) {
        jacobian = nullptr;
        Qdense = nullptr;
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
template <class AnalyticalFluxImp>
class GodunovLocalNumericalCouplingFlux
    : public LocalNumericalCouplingFluxInterface<internal::GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp>>
{
public:
  typedef internal::GodunovLocalNumericalCouplingFluxTraits<AnalyticalFluxImp> Traits;
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
template <class AnalyticalFluxImp, class BoundaryValueImp>
class GodunovLocalDirichletNumericalBoundaryFlux
    : public LocalNumericalBoundaryFluxInterface<internal::
                                                     GodunovLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp,
                                                                                                      BoundaryValueImp>>
{
public:
  typedef internal::GodunovLocalDirichletNumericalBoundaryFluxTraits<AnalyticalFluxImp, BoundaryValueImp> Traits;
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

  explicit GodunovLocalDirichletNumericalBoundaryFlux(const AnalyticalFluxType& analytical_flux,
                                                      const BoundaryValueType& boundary_values,
                                                      const XT::Common::Parameter& param,
                                                      const bool is_linear = false)
    : boundary_values_(boundary_values)
    , implementation_(analytical_flux, param, is_linear, true)
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
    const auto x_in_inside_coords = intersection.geometryInInside().global(x_in_intersection_coords);
    const RangeType u_i = local_source_entity.evaluate(x_in_inside_coords);
    const RangeType u_j = std::get<1>(local_functions_tuple)->evaluate(x_in_inside_coords);
    return implementation_.evaluate(local_functions_tuple,
                                    intersection,
                                    x_in_intersection_coords,
                                    x_in_inside_coords,
                                    x_in_inside_coords,
                                    u_i,
                                    u_j);
  } // RangeType evaluate(...) const

  // clear static variables
  static void reset()
  {
    internal::GodunovFluxImplementation<Traits>::reset();
  }

private:
  const BoundaryValueType& boundary_values_;
  const internal::GodunovFluxImplementation<Traits> implementation_;
}; // class GodunovLocalDirichletNumericalBoundaryFlux


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FLUXES_GODUNOV_HH
