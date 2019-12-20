// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_LOCAL_NUMERICAL_FLUXES_KINETIC_HH
#define DUNE_GDT_LOCAL_NUMERICAL_FLUXES_KINETIC_HH

#include <dune/xt/la/container.hh>

#include <dune/gdt/test/momentmodels/basisfunctions.hh>
#include <dune/gdt/test/momentmodels/entropyflux.hh>

#include "interface.hh"

namespace Dune {
namespace GDT {


template <class GV, class MomentBasis, class EntropyFluxImp = EntropyBasedFluxFunction<GV, MomentBasis>>
class NumericalKineticFlux
  : public NumericalFluxInterface<XT::Grid::extract_intersection_t<GV>,
                                  MomentBasis::dimFlux,
                                  MomentBasis::dimRange,
                                  typename MomentBasis::RangeFieldType>
{
  static_assert(std::is_base_of<MomentBasisInterface<typename MomentBasis::D,
                                                     MomentBasis::d,
                                                     typename MomentBasis::R,
                                                     MomentBasis::r,
                                                     MomentBasis::rC,
                                                     MomentBasis::d_flux,
                                                     MomentBasis::entropy>,
                                MomentBasis>::value,
                "Basisfunctions have to be derived from MomentBasisInterface");
  static const size_t d = MomentBasis::dimFlux;
  static const size_t m = MomentBasis::r;
  using R = typename MomentBasis::R;
  using I = XT::Grid::extract_intersection_t<GV>;
  using ThisType = NumericalKineticFlux;
  using BaseType = NumericalFluxInterface<I, d, m, R>;
  using EntropyFluxType = EntropyFluxImp;
  using SparseMatrixType = typename XT::LA::CommonSparseMatrix<R>;

public:
  using typename BaseType::DynamicStateType;
  using typename BaseType::FluxType;
  using typename BaseType::LocalIntersectionCoords;
  using typename BaseType::PhysicalDomainType;
  using typename BaseType::StateType;
  using typename BaseType::XIndependentFluxType;

  NumericalKineticFlux(const FluxType& flx, const MomentBasis& basis)
    : BaseType(flx)
    , basis_(basis)
  {}

  NumericalKineticFlux(const XIndependentFluxType& flx, const MomentBasis& basis)
    : BaseType(flx)
    , basis_(basis)
  {}

  NumericalKineticFlux(const ThisType& other) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply;
  using BaseType::flux;
  using BaseType::intersection;

  StateType apply(const LocalIntersectionCoords& /*x*/,
                  const StateType& u,
                  const StateType& v,
                  const PhysicalDomainType& n,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    StateType ret;
    apply_imp(u, v, n, ret);
    return ret;
  }

  void apply(const LocalIntersectionCoords& /*x*/,
             const DynamicStateType& u,
             const DynamicStateType& v,
             const PhysicalDomainType& n,
             DynamicStateType& ret,
             const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    apply_imp(u, v, n, ret);
  }

private:
  template <class StateTp, class RetType>
  void apply_imp(const StateTp& u, const StateTp& v, const PhysicalDomainType& n, RetType& ret) const
  {
    // find direction of unit outer normal (we assume an axis-aligned cube grid)
    size_t direction = intersection().indexInInside() / 2;
    if (dynamic_cast<const EntropyFluxType*>(&flux()) != nullptr) {
      dynamic_cast<const EntropyFluxType*>(&flux())->evaluate_kinetic_flux(
          intersection().inside(),
          intersection().neighbor() ? intersection().outside() : intersection().inside(),
          u,
          v,
          n,
          direction,
          ret);
      for (auto&& entry : ret)
        if (std::isnan(entry) || std::isinf(entry))
          DUNE_THROW(Dune::MathError, "NaN or inf in kinetic flux");
    } else {
      static const auto flux_matrices = initialize_flux_matrices(basis_);
      auto tmp_vec = ret;
      const auto& inner_flux_matrix = flux_matrices[direction][n[direction] > 0 ? 1 : 0];
      const auto& outer_flux_matrix = flux_matrices[direction][n[direction] > 0 ? 0 : 1];
      inner_flux_matrix.mv(u, tmp_vec);
      outer_flux_matrix.mv(v, ret);
      ret += tmp_vec;
      ret *= n[direction];
    }
  }

  static FieldVector<FieldVector<SparseMatrixType, 2>, d> initialize_flux_matrices(const MomentBasis& basis)
  {
    // calculate < v_i b b^T >_- M^{-1} and < v_i b b^T >_+ M^{-1}
    const auto flux_matrices_dense = basis.kinetic_flux_matrices();
    FieldVector<FieldVector<SparseMatrixType, 2>, d> flux_matrices(
        FieldVector<SparseMatrixType, 2>(SparseMatrixType(m, m, size_t(0))));
    for (size_t dd = 0; dd < d; ++dd)
      for (size_t kk = 0; kk < 2; ++kk)
        flux_matrices[dd][kk] = flux_matrices_dense[dd][kk];
    return flux_matrices;
  }

  const MomentBasis& basis_;
}; // class NumericalKineticFlux


template <class I, class MomentBasis>
NumericalKineticFlux<I, MomentBasis> make_numerical_kinetic_flux(
    const XT::Functions::
        FluxFunctionInterface<I, MomentBasis::r, MomentBasis::dimFlux, MomentBasis::r, typename MomentBasis::R>& flux,
    const MomentBasis& basis)
{
  return NumericalKineticFlux<I, MomentBasis>(flux, basis);
}

template <class I, class MomentBasis>
NumericalKineticFlux<I, MomentBasis> make_numerical_kinetic_flux(
    const XT::Functions::
        FunctionInterface<MomentBasis::r, MomentBasis::dimFlux, MomentBasis::r, typename MomentBasis::R>& flux,
    const MomentBasis& basis)
{
  return NumericalKineticFlux<I, MomentBasis>(flux, basis);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_NUMERICAL_FLUXES_KINETIC_HH
