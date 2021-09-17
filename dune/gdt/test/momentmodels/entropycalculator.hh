// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2018)

#ifndef DUNE_GDT_MOMENTMODELS_ENTROPYCALCULATOR_HH
#define DUNE_GDT_MOMENTMODELS_ENTROPYCALCULATOR_HH

#if HAVE_DUNE_XT_DATA

#  include <string>
#  include <numeric>

#  include <dune/xt/grid/functors/interfaces.hh>
#  include <dune/xt/common/parameter.hh>

#  include <dune/gdt/discretefunction/default.hh>
#  include <dune/gdt/test/momentmodels/entropyflux.hh>
#  include <dune/gdt/operators/interfaces.hh>
#  include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


template <class SpaceType, class VectorType, class MomentBasis>
class LocalEntropyCalculator : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
{
  using GridViewType = typename SpaceType::GridViewType;
  using BaseType = XT::Grid::ElementFunctor<GridViewType>;
  using EntityType = typename GridViewType::template Codim<0>::Entity;
  using IndexSetType = typename GridViewType::IndexSet;
  using EntropyFluxType = EntropyBasedFluxFunction<GridViewType, MomentBasis>;
  using RangeFieldType = typename EntropyFluxType::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;
  static constexpr size_t dimFlux = EntropyFluxType::dimFlux;
  static constexpr size_t dimRange = EntropyFluxType::basis_dimRange;
  using DiscreteFunctionType = DiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;
  using ConstDiscreteFunctionType = ConstDiscreteFunction<VectorType, GridViewType, dimRange, 1, RangeFieldType>;

public:
  explicit LocalEntropyCalculator(const SpaceType& space,
                                  const VectorType& u_n_dofs,
                                  const VectorType& u_update_dofs,
                                  const EntropyFluxType& analytical_flux,
                                  const XT::Common::Parameter& param,
                                  double& relaxationupdate)
    : space_(space)
    , u_n_(space_, u_n_dofs, "u_n")
    , u_update_(space_, u_update_dofs, "u_update")
    , local_u_n_(u_n_.local_discrete_function())
    , local_u_update_(u_update_.local_discrete_function())
    , analytical_flux_(analytical_flux)
    , local_flux_(analytical_flux_.derived_local_function())
    , param_(param)
    , index_set_(space_.grid_view().indexSet())
    , relaxationupdate_(relaxationupdate)
  {}

  explicit LocalEntropyCalculator(LocalEntropyCalculator& other)
    : BaseType(other)
    , space_(other.space_)
    , u_n_(space_, other.u_n_.dofs().vector(), "u_n")
    , u_update_(space_, other.u_update_.dofs().vector(), "u_update")
    , local_u_n_(u_n_.local_discrete_function())
    , local_u_update_(u_update_.local_discrete_function())
    , analytical_flux_(other.analytical_flux_)
    , local_flux_(analytical_flux_.derived_local_function())
    , param_(other.param_)
    , index_set_(space_.grid_view().indexSet())
    , relaxationupdate_(other.relaxationupdate_)
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalEntropyCalculator(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_u_n_->bind(entity);
    local_u_update_->bind(entity);
    XT::Common::FieldVector<RangeFieldType, dimRange> u_n;
    XT::Common::FieldVector<RangeFieldType, dimRange> u_update;
    for (size_t ii = 0; ii < dimRange; ++ii) {
      u_n[ii] = local_u_n_->dofs().get_entry(ii);
      u_update[ii] = local_u_update_->dofs().get_entry(ii);
    }
    local_flux_->bind(entity);
    const auto alpha = local_flux_->get_alpha(u_n, true)->first;
    relaxationupdate_ += std::inner_product(alpha.begin(), alpha.end(), u_update.begin(), 0.);
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType u_n_;
  const ConstDiscreteFunctionType u_update_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_u_n_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_u_update_;
  const EntropyFluxType& analytical_flux_;
  std::unique_ptr<typename EntropyFluxType::Localfunction> local_flux_;
  const XT::Common::Parameter& param_;
  const typename SpaceType::GridViewType::IndexSet& index_set_;
  double& relaxationupdate_;
}; // class LocalEntropyCalculator<...>


template <class MomentBasisImp,
          class SpaceImp,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class EntropyCalculator
  : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>
{
  using BaseType = OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>;

public:
  using typename BaseType::VectorType;
  using MomentBasis = MomentBasisImp;
  using SpaceType = SpaceImp;
  using SourceSpaceType = SpaceImp;
  using RangeSpaceType = SpaceImp;
  using EntropyFluxType = EntropyBasedFluxFunction<typename SpaceType::GridViewType, MomentBasis>;
  using RangeFieldType = typename MomentBasis::RangeFieldType;
  using LocalVectorType = typename EntropyFluxType::VectorType;

  EntropyCalculator(const EntropyFluxType& analytical_flux, const SpaceType& space)
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

  using BaseType::apply;

  void
  apply(const VectorType& /*u_n*/, VectorType& /*range*/, const XT::Common::Parameter& /*param*/) const override final
  {
    DUNE_THROW(Dune::NotImplemented, "Kaputt!");
  } // void apply(...)

  void
  apply(const VectorType& u_n, VectorType& u_update, const XT::Common::Parameter& param, double& relaxationupdate) const
  {
    LocalEntropyCalculator<SpaceType, VectorType, MomentBasis> local_entropy_calculator(
        space_, u_n, u_update, analytical_flux_, param, relaxationupdate);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_entropy_calculator);
    walker.walk(false);
  } // void apply(...)

private:
  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
}; // class EntropyCalculator<...>


} // namespace GDT
} // namespace Dune

#endif // HAVE_DUNE_XT_DATA

#endif // DUNE_GDT_MOMENTMODELS_ENTROPYCALCULATOR_HH
