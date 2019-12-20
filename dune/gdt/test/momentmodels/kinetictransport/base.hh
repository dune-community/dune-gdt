// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
#define DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH

#include <dune/grid/common/partitionset.hh>

#include <dune/xt/functions/constant.hh>

#include <dune/xt/la/solver.hh>

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>
#include <dune/gdt/test/momentmodels/kineticequation.hh>

#include "../kineticequation.hh"

namespace Dune {
namespace GDT {


template <class E, class MomentBasisImp>
class KineticTransportEquationBase : public KineticEquationInterface<E, MomentBasisImp>
{
  using ThisType = KineticTransportEquationBase;
  using BaseType = KineticEquationInterface<E, MomentBasisImp>;

public:
  using BaseType::dimDomain;
  using BaseType::dimFlux;
  using BaseType::dimRange;
  using BaseType::dimRangeCols;
  using typename BaseType::BasisDomainType;
  using typename BaseType::DomainFieldType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::GenericFluxFunctionType;
  using typename BaseType::GenericFunctionType;
  using typename BaseType::MatrixType;
  using typename BaseType::MomentBasis;
  using typename BaseType::RangeFieldType;
  using typename BaseType::RangeReturnType;
  using typename BaseType::StateType;

  using typename BaseType::BoundaryValueType;
  using typename BaseType::FluxType;
  using typename BaseType::InitialValueType;
  using FluxRangeType = typename FluxType::LocalFunctionType::RangeReturnType;
  using DynamicFluxRangeType = typename FluxType::LocalFunctionType::DynamicRangeType;
  using FluxJacobianRangeType = typename FluxType::LocalFunctionType::JacobianRangeReturnType;
  using DynamicFluxJacobianRangeType = typename FluxType::LocalFunctionType::DynamicJacobianRangeType;
  using GenericScalarFunctionType = XT::Functions::GenericFunction<dimFlux, 1, 1, RangeFieldType>;
  using ConstantScalarFunctionType = XT::Functions::ConstantFunction<dimFlux, 1, 1, RangeFieldType>;
  using BoundaryDistributionType =
      std::function<std::function<RangeFieldType(const BasisDomainType&)>(const DomainType&)>;

  using BaseType::default_boundary_cfg;
  using BaseType::default_grid_cfg;

  KineticTransportEquationBase(const MomentBasis& basis_functions,
                               const XT::Common::Configuration& grid_cfg = default_grid_cfg(),
                               const XT::Common::Configuration& boundary_cfg = default_boundary_cfg(),
                               const RangeFieldType psi_vac = 5e-9)
    : BaseType(basis_functions)
    , grid_cfg_(grid_cfg)
    , boundary_cfg_(boundary_cfg)
    , psi_vac_(psi_vac)
  {}

  template <class VectorType>
  void solve(const MatrixType& mat,
             VectorType& x,
             const VectorType& rhs,
             const MomentBasisInterface<DomainFieldType,
                                        MomentBasis::dimDomain,
                                        RangeFieldType,
                                        dimRange,
                                        dimRangeCols,
                                        dimFlux,
                                        MomentBasis::entropy>&) const
  {
    // copy to CommonDenseMatrix as the FieldMatrix copies itself on the stack during solve which may case a
    // stackoverflow for large matrices
    XT::LA::CommonDenseMatrix<RangeFieldType> xt_la_mat(mat);
    XT::LA::CommonDenseVector<RangeFieldType> xt_la_rhs(rhs, 0);
    XT::LA::CommonDenseVector<RangeFieldType> xt_la_x(rhs.size());
    XT::LA::solve(xt_la_mat, xt_la_rhs, xt_la_x);
    for (size_t ii = 0; ii < xt_la_x.size(); ++ii)
      x[ii] = xt_la_x[ii];
  }

  template <class VectorType, size_t refinements>
  void solve(const MatrixType& mat,
             VectorType& x,
             const VectorType& rhs,
             const PartialMomentBasis<DomainFieldType, 3, RangeFieldType, refinements, dimRangeCols, 3, 1>&) const
  {
    const size_t num_blocks = dimRange / 4;
    const size_t block_size = 4;
    XT::Common::FieldMatrix<RangeFieldType, block_size, block_size> local_mat;
    XT::Common::FieldVector<RangeFieldType, block_size> local_x, local_rhs;
    for (size_t jj = 0; jj < num_blocks; ++jj) {
      const size_t offset = jj * block_size;
      // copy to local matrix and vector
      for (size_t rr = 0; rr < block_size; ++rr) {
        local_rhs[rr] = rhs[offset + rr];
        for (size_t cc = 0; cc < block_size; ++cc)
          local_mat[rr][cc] = mat[offset + rr][offset + cc];
      } // rr
      local_mat.solve(local_x, local_rhs);
      for (size_t rr = 0; rr < block_size; ++rr)
        x[offset + rr] = local_x[rr];
    } // jj
  }

  // flux matrix A = B M^{-1} with B_{ij} = <v h_i h_j>
  std::unique_ptr<FluxType> flux() const override
  {
    // calculate B row-wise by solving M^{T} A^T = B^T column-wise
    auto A = basis_functions_.flux_matrix();
    const auto M_T = basis_functions_.mass_matrix(); // mass matrix is symmetric
    // solve
    DynamicVector<RangeFieldType> tmp_row(M_T.N(), 0.);
    for (size_t dd = 0; dd < dimFlux; ++dd) {
      for (size_t ii = 0; ii < M_T.N(); ++ii) {
        solve(M_T, tmp_row, A[dd][ii], basis_functions_);
        A[dd][ii] = tmp_row;
      }
    }
    auto order_func = [](const XT::Common::Parameter&) -> int { return 1; };
    DynamicVector<XT::LA::CommonDenseMatrix<RangeFieldType>> A_la(dimFlux);
    for (size_t dd = 0; dd < dimFlux; ++dd)
      A_la[dd] = A[dd];
    auto eval_func =
        [A_la](const DomainType&, const StateType& u, DynamicFluxRangeType& ret, const XT::Common::Parameter&) {
          for (size_t dd = 0; dd < dimFlux; ++dd) {
            auto row_view = ret[dd];
            A_la[dd].mv(u, row_view);
          }
        };
    auto jacobian_func =
        [A_la](const DomainType&, const StateType&, DynamicFluxJacobianRangeType& ret, const XT::Common::Parameter&) {
          for (size_t dd = 0; dd < dimFlux; ++dd)
            ret[dd] = A_la[dd];
        };
    return std::make_unique<GenericFluxFunctionType>(order_func,
                                                     GenericFluxFunctionType::default_post_bind_function(),
                                                     eval_func,
                                                     XT::Common::ParameterType{},
                                                     "flux",
                                                     jacobian_func);
  }

  // Initial value of the kinetic equation is a constant vacuum concentration psi_vac.
  // Thus, the initial value of the n-th moment is basis_integrated * psi_vac.
  std::unique_ptr<InitialValueType> initial_values() const override
  {
    const auto value = basis_functions_.integrated() * psi_vac_;
    return std::make_unique<GenericFunctionType>(
        [](const XT::Common::Parameter&) { return 0; },
        [value](const DomainType&, DynamicRangeType& ret, const XT::Common::Parameter&) { ret = value; });
  } // ... initial_values()

  // Use a constant vacuum concentration basis_integrated * psi_vac as default boundary value
  std::unique_ptr<BoundaryValueType> boundary_values() const override
  {
    const auto value = basis_functions_.integrated() * psi_vac_;
    return std::make_unique<GenericFunctionType>(
        [](const XT::Common::Parameter&) { return 0; },
        [=](const DomainType&, DynamicRangeType& ret, const XT::Common::Parameter&) { ret = value; });
  } // ... boundary_values()

  virtual BoundaryDistributionType boundary_distribution() const
  {
    return [this](const DomainType&) { return [this](const BasisDomainType&) { return this->psi_vac_; }; };
  }

  RangeReturnType
  kinetic_boundary_flux_from_quadrature(const DomainType& x, const RangeFieldType& n, const size_t dd) const
  {
    RangeReturnType ret(0.);
    const auto boundary_density = boundary_distribution()(x);
    const auto& quadratures = basis_functions_.quadratures();
    for (size_t jj = 0; jj < quadratures.size(); ++jj) {
      for (size_t ll = 0; ll < quadratures[jj].size(); ++ll) {
        const auto v = quadratures[jj][ll].position();
        const auto b = basis_functions_.evaluate(v, jj);
        if (v[dd] * n < 0) {
          const RangeFieldType psi = boundary_density(v);
          ret += b * psi * v[dd] * quadratures[jj][ll].weight();
        }
      } // ll
    } // jj
    return ret;
  }

  virtual RangeReturnType kinetic_boundary_flux(const DomainType& x, const RangeFieldType& n, const size_t dd) const
  {
    return kinetic_boundary_flux_from_quadrature(x, n, dd);
  }

  RangeFieldType CFL() const override
  {
    return 0.49;
  }

  RangeFieldType t_end() const override
  {
    return 1.;
  }

  XT::Common::Configuration grid_config() const override
  {
    return grid_cfg_;
  }

  XT::Common::Configuration boundary_config() const override
  {
    return boundary_cfg_;
  }

  static std::string static_id()
  {
    return "kinetictransportequation";
  }

  virtual const RangeFieldType psi_vac() const
  {
    return psi_vac_;
  }

protected:
  using BaseType::basis_functions_;
  const XT::Common::Configuration grid_cfg_;
  const XT::Common::Configuration boundary_cfg_;
  const RangeFieldType psi_vac_;
  XT::Common::ParameterType parameter_type_;
}; // class KineticTransportEquation<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_HYPERBOLIC_PROBLEMS_MOMENTMODELS_KINETICTRANSPORTEQUATION_BASE_HH
