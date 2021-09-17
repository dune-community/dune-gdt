// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Rene Milk      (2017 - 2018)
//   Tobias Leibner (2017)

#ifndef DUNE_GDT_OPERATORS_FV_RHS_HH
#define DUNE_GDT_OPERATORS_FV_RHS_HH

#include <boost/config.hpp>

#include <dune/xt/common/mkl.hh>
#include <dune/xt/common/numeric.hh>
#include <dune/gdt/operators/interfaces.hh>

namespace Dune {
namespace GDT {

template <class SpaceType, class VectorType, class MomentBasis>
class LocalRhsOperator : public XT::Grid::ElementFunctor<typename SpaceType::GridViewType>
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
  using RangeFieldVector = XT::Common::FieldVector<RangeFieldType, dimRange>;
  using ScalarFunctionType = XT::Functions::FunctionInterface<dimFlux, 1, 1, RangeFieldType>;

public:
  explicit LocalRhsOperator(const SpaceType& space,
                            const VectorType& source_dofs,
                            VectorType& range_dofs,
                            const MomentBasis& basis_functions,
                            const RangeFieldVector& basis_integrated,
                            const RangeFieldVector& u_iso,
                            const ScalarFunctionType& sigma_a,
                            const ScalarFunctionType& sigma_s,
                            const ScalarFunctionType& Q)
    : space_(space)
    , source_(space_, source_dofs, "source")
    , local_source_(source_.local_discrete_function())
    , range_(space_, range_dofs, "range")
    , local_range_(range_.local_discrete_function())
    , basis_functions_(basis_functions)
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
    , u_iso_(u_iso)
    , basis_integrated_(basis_integrated)
    , index_set_(space_.grid_view().indexSet())
  {}

  explicit LocalRhsOperator(LocalRhsOperator& other)
    : BaseType(other)
    , space_(other.space_)
    , source_(space_, other.source_.dofs().vector(), "source")
    , local_source_(source_.local_discrete_function())
    , range_(space_, other.range_.dofs().vector(), "range")
    , local_range_(range_.local_discrete_function())
    , basis_functions_(other.basis_functions_)
    , sigma_a_(other.sigma_a_)
    , sigma_s_(other.sigma_s_)
    , Q_(other.Q_)
    , u_iso_(other.u_iso_)
    , basis_integrated_(other.basis_integrated_)
    , index_set_(space_.grid_view().indexSet())
  {}

  XT::Grid::ElementFunctor<GridViewType>* copy() override final
  {
    return new LocalRhsOperator(*this);
  }

  void apply_local(const EntityType& entity) override final
  {
    local_source_->bind(entity);
    local_range_->bind(entity);
    XT::Common::FieldVector<RangeFieldType, dimRange> u;
    for (size_t ii = 0; ii < dimRange; ++ii)
      u[ii] = local_source_->dofs().get_entry(ii);
    const auto center = entity.geometry().center();
    const auto sigma_a_value = sigma_a_.evaluate(center)[0];
    const auto sigma_s_value = sigma_s_.evaluate(center)[0];
    const auto sigma_t_value = sigma_a_value + sigma_s_value;
    const auto Q_value = Q_.evaluate(center)[0];
    auto ret = u;
    ret *= -sigma_t_value;
    ret.axpy(basis_functions_.density(u) * sigma_s_value, u_iso_);
    ret.axpy(Q_value, basis_integrated_);
    for (size_t ii = 0; ii < dimRange; ++ii)
      local_range_->dofs().add_to_entry(ii, ret[ii]);
  } // void apply_local(...)

private:
  const SpaceType& space_;
  const ConstDiscreteFunctionType source_;
  std::unique_ptr<typename ConstDiscreteFunctionType::ConstLocalDiscreteFunctionType> local_source_;
  DiscreteFunctionType range_;
  std::unique_ptr<typename DiscreteFunctionType::LocalDiscreteFunctionType> local_range_;
  const MomentBasis& basis_functions_;
  const ScalarFunctionType& sigma_a_;
  const ScalarFunctionType& sigma_s_;
  const ScalarFunctionType& Q_;
  const RangeFieldVector& u_iso_;
  const RangeFieldVector& basis_integrated_;
  const typename SpaceType::GridViewType::IndexSet& index_set_;
}; // class LocalRhsOperator<...>

template <class MomentBasisImp,
          class SpaceImp,
          class MatrixType = typename XT::LA::Container<typename MomentBasisImp::RangeFieldType>::MatrixType>
class RhsOperator : public OperatorInterface<MatrixType, typename SpaceImp::GridViewType, MomentBasisImp::dimRange, 1>
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
  static constexpr size_t dimFlux = EntropyFluxType::dimFlux;
  static constexpr size_t dimRange = EntropyFluxType::basis_dimRange;
  using ScalarFunctionType = XT::Functions::FunctionInterface<dimFlux, 1, 1, RangeFieldType>;
  using RangeFieldVector = XT::Common::FieldVector<RangeFieldType, dimRange>;

  RhsOperator(const EntropyFluxType& analytical_flux,
              const SpaceType& space,
              const ScalarFunctionType& sigma_a,
              const ScalarFunctionType& sigma_s,
              const ScalarFunctionType& Q)
    : analytical_flux_(analytical_flux)
    , space_(space)
    , basis_functions_(analytical_flux_.basis_functions())
    , sigma_a_(sigma_a)
    , sigma_s_(sigma_s)
    , Q_(Q)
    , u_iso_(basis_functions_.u_iso())
    , basis_integrated_(basis_functions_.integrated())
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

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/) const override final
  {
    LocalRhsOperator<SpaceType, VectorType, MomentBasis> local_rhs_operator(
        space_, source, range, basis_functions_, basis_integrated_, u_iso_, sigma_a_, sigma_s_, Q_);
    auto walker = XT::Grid::Walker<typename SpaceType::GridViewType>(space_.grid_view());
    walker.append(local_rhs_operator);
    walker.walk(true);
  }

  const EntropyFluxType& analytical_flux_;
  const SpaceType& space_;
  const MomentBasis& basis_functions_;
  const ScalarFunctionType& sigma_a_;
  const ScalarFunctionType& sigma_s_;
  const ScalarFunctionType& Q_;
  const RangeFieldVector u_iso_;
  const RangeFieldVector basis_integrated_;
}; // class RhsOperator<...>


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_FV_RHS_HH
