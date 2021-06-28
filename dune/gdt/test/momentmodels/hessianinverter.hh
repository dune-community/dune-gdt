// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_MOMENTMODELS_HESSIANINVERTER_HH
#define DUNE_GDT_MOMENTMODELS_HESSIANINVERTER_HH

#if HAVE_DUNE_XT_DATA

#  include <string>

#  include <dune/xt/grid/functors/interfaces.hh>
#  include <dune/xt/common/parameter.hh>

#  include <dune/gdt/discretefunction/default.hh>
#  include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#  include <dune/gdt/operators/interfaces.hh>
#  include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class MomentBasis, SlopeLimiterType slope>
class LocalEntropicHessianInverter : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using GridViewType = typename SpaceType::GridViewType;
  using BaseType = XT::Grid::ElementFunctor<GridViewType>;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedFluxEntropyCoordsFunction<GridViewType, MomentBasis, slope>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  static constexpr size_t dimFlux = EntropyFluxType::dimFlux;
  static constexpr size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;

public:
  explicit LocalEntropicHessianInverter(const SpaceType& space,
                                        const VectorType& u_update_dofs,
                                        const VectorType& alpha_dofs,
                                        std::vector<bool>& reg_indicators,
                                        VectorType& alpha_range_dofs,
                                        const EntropyFluxType& analytical_flux,
                                        const XT::Common::Parameter& param,
                                        double& relaxationupdate)
    : space_(space)
    , u_update_(space_, u_update_dofs, "u_update")
    , alpha_(space_, alpha_dofs, "u_update")
    , reg_indicators_(reg_indicators)
    , local_u_update_(u_update_.local_discrete_function())
    , local_alpha_(alpha_.local_discrete_function())
    , range_(space_, alpha_range_dofs, "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , param_(param)
    , relaxationupdate_(relaxationupdate)
  {}

  explicit LocalEntropicHessianInverter(LocalEntropicHessianInverter& other)
    : BaseType(other)
    , space_(other.space_)
    , u_update_(space_, other.u_update_.dofs().vector(), "source")
    , alpha_(space_, other.alpha_.dofs().vector(), "source")
    , reg_indicators_(other.reg_indicators_)
    , local_u_update_(u_update_.local_discrete_function())
    , local_alpha_(alpha_.local_discrete_function())
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , param_(other.param_)
    , relaxationupdate_(other.relaxationupdate_)
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalEntropicHessianInverter(*this);
  }

  /** \short applies inverse Hessian
   *
   *  Returns H^{-1} u, where H is the Hessian on this entity and u is the update vector in original coordinates
   *  If inversion fails, an error is thrown. If regularization is requested (by setting the key "reg" in param_), the
   *  entity is marked for regularization instead.
   *  \note Evaluations of the ansatz densities have to be stored in
   *  analytical_flux_ (by calling store_evaluations on the analytical flux before using this function, see
   *  density_evaluator.hh).
   *  \note Regularization is poorly tested and probably does not really work, a lot of tuning is
   *  still missing (when do we want to use regularization, and when is it better to simply reduce the timestep?).
   */
  void apply_local(const EntityType& entity) override final
  {
    local_u_update_->bind(entity);
    local_alpha_->bind(entity);
    local_range_->bind(entity);
    const auto& local_u_dofs = local_u_update_->dofs();
    const auto& local_alpha_dofs = local_alpha_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      Hinv_u_[ii] = local_u_dofs.get_entry(ii);
    for (size_t ii = 0; ii < dimRange; ++ii)
      relaxationupdate_ += local_alpha_dofs.get_entry(ii) * Hinv_u_[ii];
    const auto entity_index = space_.grid_view().indexSet().index(entity);
    const double dt = param_.get("dt")[0];
    try {
      analytical_flux_.apply_inverse_hessian(entity_index, Hinv_u_, dt);
      for (auto&& entry : Hinv_u_)
        if (std::isnan(entry) || std::isinf(entry))
          DUNE_THROW(Dune::MathError, "Hessian");
    } catch (const Dune::MathError& e) {
      if (param_.has_key("reg") && param_.get("reg")[0]) {
        std::cout << "reg considered" << std::endl;
        for (size_t ii = 0; ii < dimRange; ++ii)
          Hinv_u_[ii] = local_u_dofs.get_entry(ii);
        const auto rho = analytical_flux_.basis_functions().density(analytical_flux_.get_u(Hinv_u_));
        if ((rho < 1e-7 && dt < 1e-4) || dt < 1e-7) {
          std::cout << "reg indicator set" << std::endl;
          reg_indicators_[entity_index] = true;
        }
        return;
      } else
        throw e;
    }
    auto& local_range_dofs = local_range_->dofs();
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_dofs.set_entry(ii, Hinv_u_[ii]);
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType u_update_;
  const ConstDiscreteFunctionType alpha_;
  std::vector<bool>& reg_indicators_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_u_update_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_alpha_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  const EntropyFluxType& analytical_flux_;
  const XT::Common::Parameter& param_;
  XT::Common::FieldVector<RangeFieldType, dimRange> u_tmp_;
  XT::Common::FieldVector<RangeFieldType, dimRange> Hinv_u_;
  double& relaxationupdate_;
}; // class LocalEntropicHessianInverter<...>

template <class MomentBasisImp,
          class SpaceImp,
          SlopeLimiterType slope,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class EntropicHessianInverter
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

  EntropicHessianInverter(const EntropyFluxType& analytical_flux, const SpaceType& space)
    : analytical_flux_(analytical_flux)
    , space_(space)
  {}

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

  void apply(const VectorType& /*source*/,
             VectorType& /*range*/,
             const XT::Common::Parameter& /*param*/) const override final
  {
    DUNE_THROW(Dune::NotImplemented, "Use apply_inverse_hessian!");
  } // void apply(...)

  void apply_inverse_hessian(const VectorType& u_update,
                             const VectorType& alpha,
                             std::vector<bool>& reg_indicators,
                             VectorType& alpha_update,
                             const XT::Common::Parameter& param,
                             double& relaxationupdate) const
  {
    LocalEntropicHessianInverter<SpaceType, VectorType, MomentBasis, slope> local_hessian_inverter(
        space_, u_update, alpha, reg_indicators, alpha_update, analytical_flux_, param, relaxationupdate);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_hessian_inverter);
    walker.walk(true);
  } // void apply(...)

  //  template <class ElementRange>
  //  void apply_inverse_hessian_range(const VectorType& u_update,
  //                                   std::vector<bool>& reg_indicators,
  //                                   VectorType& alpha_update,
  //                                   const XT::Common::Parameter& param,
  //                                   const ElementRange& element_range) const
  //  {
  //    LocalEntropicHessianInverter<SpaceType, VectorType, MomentBasis, slope> local_hessian_inverter(
  //        space_, u_update, alpha, reg_indicators, alpha_update, analytical_flux_, param);
  //    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
  //    walker.append(local_hessian_inverter);
  //    walker.walk_range(element_range);
  //  } // void apply(...)


private:
  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
}; // class EntropicHessianInverter<...>


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_MOMENTMODELS_HESSIANINVERTER_HH
