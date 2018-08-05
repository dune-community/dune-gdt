// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_ADVECTION_FV_HH
#define DUNE_GDT_OPERATORS_ADVECTION_FV_HH

#include <dune/xt/common/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/operators/advection-fv.hh>

#include "interfaces.hh"
#include "localizable-operator.hh"

namespace Dune {
namespace GDT {


/**
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class AssemblyGridView,
          class SV,
          size_t m = 1,
          class SF = double,
          class SGV = AssemblyGridView,
          class RF = SF,
          class RGV = AssemblyGridView,
          class RV = SV>
class AdvectionFvOperator : public OperatorInterface<SV,
                                                     SGV,
                                                     m,
                                                     1,
                                                     SF,
                                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                                     m,
                                                     1,
                                                     RF,
                                                     RGV,
                                                     RV>
{
  // No need to check the rest, is done in OperatorInterface.
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert((AssemblyGridView::dimension == SGV::dimension) && (AssemblyGridView::dimension == RGV::dimension), "");

  using ThisType = AdvectionFvOperator<AssemblyGridView, SV, m, SF, SGV, RF, RGV, RV>;
  using BaseType = OperatorInterface<SV,
                                     SGV,
                                     m,
                                     1,
                                     SF,
                                     typename XT::Common::multiplication_promotion<SF, RF>::type,
                                     m,
                                     1,
                                     RF,
                                     RGV,
                                     RV>;
  using AGV = AssemblyGridView;

public:
  using AssemblyGridViewType = AssemblyGridView;
  using NumericalFluxType = NumericalFluxInterface<AssemblyGridView::dimension, m, RF>;

  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;

  using typename BaseType::SourceVectorType;
  using typename BaseType::RangeVectorType;

  AdvectionFvOperator(const AssemblyGridViewType& assembly_grid_view,
                      const NumericalFluxType& numerical_flux,
                      const SourceSpaceType& source_space,
                      const RangeSpaceType& range_space,
                      const XT::Grid::IntersectionFilter<AssemblyGridViewType>& periodicity_exception =
                          XT::Grid::ApplyOn::NoIntersections<AssemblyGridViewType>())
    : BaseType(numerical_flux.parameter_type())
    , assembly_grid_view_(assembly_grid_view)
    , numerical_flux_(numerical_flux)
    , source_space_(source_space)
    , range_space_(range_space)
    , periodicity_exception_(periodicity_exception.copy())
  {
  }

  AdvectionFvOperator(ThisType&& source) = default;

  bool linear() const override final
  {
    return numerical_flux_.linear();
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  using BaseType::apply;

  void apply(const SourceVectorType& source,
             RangeVectorType& range,
             const XT::Common::Parameter& param = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
    const auto source_function = make_discrete_function(source_space_, source);
    auto range_function = make_discrete_function(range_space_, range);
    // set up the actual operator
    auto localizable_op = make_localizable_operator(assembly_grid_view_, source_function, range_function);
    using I = XT::Grid::extract_intersection_t<AGV>;
    // treat all inner intersections
    localizable_op.append(LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(numerical_flux_),
                          param,
                          XT::Grid::ApplyOn::InnerIntersectionsOnce<AGV>());
    // treat periodic boundaries
    localizable_op.append(LocalAdvectionFvCouplingOperator<I, SV, SGV, m, SF, RF, RGV, RV>(numerical_flux_),
                          param,
                          *(XT::Grid::ApplyOn::PeriodicBoundaryIntersectionsOnce<AGV>() && !(*periodicity_exception_)));
    // do the actual work
    localizable_op.assemble(/*use_tbb=*/true);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

private:
  const AssemblyGridViewType& assembly_grid_view_;
  const NumericalFluxType& numerical_flux_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::unique_ptr<XT::Grid::IntersectionFilter<AssemblyGridViewType>> periodicity_exception_;
}; // class AdvectionFvOperator


template <class VectorType, class AGV, size_t m, class RF, class SGV, class SF, class RGV>
std::enable_if_t<XT::LA::is_vector<VectorType>::value,
                 AdvectionFvOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>>
make_advection_fv_operator(
    const AGV& assembly_grid_view,
    const NumericalFluxInterface<AGV::dimension, m, RF>& numerical_flux,
    const SpaceInterface<SGV, m, 1, SF>& source_space,
    const SpaceInterface<RGV, m, 1, RF>& range_space,
    const XT::Grid::IntersectionFilter<AGV>& periodicity_exception = XT::Grid::ApplyOn::NoIntersections<AGV>())
{
  return AdvectionFvOperator<AGV, VectorType, m, SF, SGV, RF, RGV, VectorType>(
      assembly_grid_view, numerical_flux, source_space, range_space, periodicity_exception);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ADVECTION_FV_HH
