// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_MOMENTMODELS_MINDENSITYSETTER_HH
#define DUNE_GDT_MOMENTMODELS_MINDENSITYSETTER_HH

#if HAVE_DUNE_XT_DATA

#  include <atomic>

#  include <dune/xt/common/parameter.hh>
#  include <dune/xt/grid/functors/interfaces.hh>

#  include <dune/gdt/discretefunction/default.hh>
#  include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#  include <dune/gdt/operators/interfaces.hh>
#  include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class MomentBasis, SlopeLimiterType slope>
class LocalMinDensitySetter : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using BaseType = XT::Grid::ElementFunctor<typename SpaceType::GridViewType>;

public:
  using GridViewType = typename SpaceType::GridViewType;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GridViewType, MomentBasis, slope>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  static constexpr size_t dimFlux = EntropyFluxType::dimFlux;
  static constexpr size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using DomainType = FieldVector<RangeFieldType, dimFlux>;

  explicit LocalMinDensitySetter(const SpaceType& space,
                                 const VectorType& alpha_dofs,
                                 VectorType& range_dofs,
                                 EntropyFluxType& analytical_flux,
                                 const RangeFieldType min_acceptable_density,
                                 const double dt,
                                 std::atomic<bool>& changed)
    : space_(space)
    , alpha_(space_, alpha_dofs, "alpha")
    , range_(space_, range_dofs, "regularized alpha")
    , local_alpha_(alpha_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , min_acceptable_density_(min_acceptable_density)
    , dt_(dt)
    , changed_(changed)
  {}

  explicit LocalMinDensitySetter(LocalMinDensitySetter& other)
    : BaseType(other)
    , space_(other.space_)
    , alpha_(space_, other.alpha_.dofs().vector(), "source")
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_alpha_(alpha_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , min_acceptable_density_(other.min_acceptable_density_)
    , dt_(other.dt_)
    , changed_(other.changed_)
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalMinDensitySetter(*this);
  }

  /** \short Modifies alpha to ensure that the density is larger than min_acceptable_density (or to improve stability,
   * in some cases)
   *
   * The actual modification is dependent on the basis function, so consider the adjust_alpha_to_ensure_min_density
   * method for the respective basis function class. Currently only really useful for the Hatfunction basis. Needs some
   * more investigation. Unused by default at the moment, to enable add -adjust_alpha 1 to the command line.
   * Modifications are only applied if the current timestep is lesser than adjust_dt, which is 1e-3 by default and can
   * be modified on the command line via -adjust_dt dt (e.g., -adjust_dt 1e-4)
   */
  void apply_local(const EntityType& entity) override final
  {
    local_alpha_->bind(entity);
    local_range_->bind(entity);
    const auto& basis_functions = analytical_flux_.basis_functions();
    const auto& local_alpha_dofs = local_alpha_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      alpha_tmp_[ii] = local_alpha_dofs.get_entry(ii);
    static const double adjust_dt = DXTC_CONFIG_GET("adjust_dt", 1.0e-3);
    if (dt_ < adjust_dt) {
      const bool changed = basis_functions.adjust_alpha_to_ensure_min_density(
          alpha_tmp_, min_acceptable_density_, analytical_flux_.get_u(alpha_tmp_));
      if (changed)
        changed_ = true;
    }
    auto& local_range_dofs = local_range_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_dofs.set_entry(ii, alpha_tmp_[ii]);
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType alpha_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_alpha_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  EntropyFluxType& analytical_flux_;
  const RangeFieldType min_acceptable_density_;
  const RangeFieldType dt_;
  XT::Common::FieldVector<RangeFieldType, dimRange> alpha_tmp_;
  std::atomic<bool>& changed_;
}; // class LocalMinDensitySetter<...>

template <class MomentBasisImp,
          class SpaceImp,
          SlopeLimiterType slope,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class MinDensitySetter
  : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>
{
  using BaseType = OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>;

public:
  using typename BaseType::VectorType;
  using MomentBasis = MomentBasisImp;
  using SpaceType = SpaceImp;
  using SourceSpaceType = SpaceImp;
  using RangeSpaceType = SpaceImp;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<typename SpaceType::GridViewType, MomentBasis, slope>;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using LocalMinDensitySetterType = LocalMinDensitySetter<SpaceType, VectorType, MomentBasis, slope>;
  using EntityType = typename LocalMinDensitySetterType::EntityType;

  MinDensitySetter(EntropyFluxType& analytical_flux,
                   const SpaceType& space,
                   const RangeFieldType min_acceptable_density)
    : analytical_flux_(analytical_flux)
    , space_(space)
    , min_acceptable_density_(min_acceptable_density)
  {
    analytical_flux_.prepare_storage(space_.grid_view());
  }

  bool linear() const override final
  {
    return false;
  }

  const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  void apply(const VectorType& alpha, VectorType& range, const XT::Common::Parameter& param = {}) const override final
  {
    static const bool adjust = DXTC_CONFIG_GET("adjust_alpha", 0);
    if (!adjust)
      return;
    std::atomic<bool> changed(false);
    LocalMinDensitySetterType local_min_density_setter(
        space_, alpha, range, analytical_flux_, min_acceptable_density_, param.get("dt")[0], changed);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_min_density_setter);
    walker.walk(true);
  } // void apply(...)

  bool apply_with_dt(const VectorType& alpha, VectorType& range, const double dt) const
  {
    static const bool adjust = DXTC_CONFIG_GET("adjust_alpha", 0);
    if (!adjust)
      return false;
    std::atomic<bool> changed(false);
    LocalMinDensitySetterType local_min_density_setter(
        space_, alpha, range, analytical_flux_, min_acceptable_density_, dt, changed);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_min_density_setter);
    walker.walk(true);
    return changed;
  } // void apply(...)

private:
  EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
  const RangeFieldType min_acceptable_density_;
}; // class MinDensitySetter<...>


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_MOMENTMODELS_MINDENSITYSETTER_HH
