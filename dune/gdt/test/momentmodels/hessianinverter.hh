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

#include <string>

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/common/parameter.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/test/momentmodels/entropyflux_kineticcoords.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/type_traits.hh>

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
  static const size_t dimFlux = EntropyFluxType::dimFlux;
  static const size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;

public:
  explicit LocalEntropicHessianInverter(const SpaceType& space,
                                        const VectorType& alpha_dofs,
                                        const VectorType& u_update_dofs,
                                        std::vector<bool>& reg_indicators,
                                        VectorType& alpha_range_dofs,
                                        const EntropyFluxType& analytical_flux,
                                        const XT::Common::Parameter& param)
    : space_(space)
    , alpha_(space_, alpha_dofs, "alpha")
    , u_update_(space_, u_update_dofs, "u_update")
    , reg_indicators_(reg_indicators)
    , local_alpha_(alpha_.local_discrete_function())
    , local_u_update_(u_update_.local_discrete_function())
    , range_(space_, alpha_range_dofs, "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , param_(param)
  {}

  explicit LocalEntropicHessianInverter(LocalEntropicHessianInverter& other)
    : BaseType(other)
    , space_(other.space_)
    , alpha_(space_, other.alpha_.dofs().vector(), "source")
    , u_update_(space_, other.u_update_.dofs().vector(), "source")
    , reg_indicators_(other.reg_indicators_)
    , local_alpha_(alpha_.local_discrete_function())
    , local_u_update_(u_update_.local_discrete_function())
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_range_(range_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , param_(other.param_)
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalEntropicHessianInverter(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_u_update_->bind(entity);
    local_range_->bind(entity);
    XT::Common::FieldVector<RangeFieldType, dimRange> u, Hinv_u;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u[ii] = local_u_update_->dofs().get_entry(ii);
    const auto entity_index = space_.grid_view().indexSet().index(entity);
    try {
      analytical_flux_.apply_inverse_hessian(entity_index, u, Hinv_u);
    } catch (const Dune::MathError& e) {
      if (param_.has_key("reg") && param_.get("reg")[0]) {
        reg_indicators_[entity_index] = true;
        return;
      } else
        throw e;
    }
    for (auto&& entry : Hinv_u)
      if (std::isnan(entry) || std::isinf(entry)) {
        //        std::cout << "x: " << entity.geometry().center() << "u: " << u << ", alpha: " << alpha << ", Hinv_u: "
        //        << Hinv_u << std::endl;
        DUNE_THROW(Dune::MathError, "Hessian");
      }

    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_->dofs().set_entry(ii, Hinv_u[ii]);
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType alpha_;
  const ConstDiscreteFunctionType u_update_;
  std::vector<bool>& reg_indicators_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_alpha_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_u_update_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  const EntropyFluxType& analytical_flux_;
  const XT::Common::Parameter& param_;
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

  void apply_inverse_hessian(const VectorType& alpha,
                             const VectorType& u_update,
                             std::vector<bool>& reg_indicators,
                             VectorType& alpha_update,
                             const XT::Common::Parameter& param) const
  {
    LocalEntropicHessianInverter<SpaceType, VectorType, MomentBasis, slope> local_hessian_inverter(
        space_, alpha, u_update, reg_indicators, alpha_update, analytical_flux_, param);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_hessian_inverter);
    walker.walk(true);
  } // void apply(...)

private:
  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
}; // class EntropicHessianInverter<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_MOMENTMODELS_HESSIANINVERTER_HH
