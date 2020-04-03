// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
#define DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH

#include <list>

#include <dune/xt/common/deprecated.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/grid/walker.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/assembler/operator-applicators.hh>
#include <dune/gdt/local/operators/generic.hh>
#include <dune/gdt/local/operators/interfaces.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \todo Create LocalizableOperatorApplicator which accepts a GridFunction as source, derive
 *       LocalizableDiscreteOperatorApplicator from LocalizableOperatorApplicator.
 *
 * \note Most likely, you want to use LocalizableOperator.
 */
template <class AssemblyGridView,
          class SourceVector,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          class SourceGridView = AssemblyGridView,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
class LocalizableDiscreteOperatorApplicator : public XT::Grid::Walker<AssemblyGridView>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_vector<SourceVector>::value, "");
  static_assert(XT::Grid::is_view<SourceGridView>::value, "");
  static_assert(XT::Grid::is_view<RangeGridView>::value, "");
  static_assert(XT::LA::is_vector<RangeVector>::value, "");

  using ThisType = LocalizableDiscreteOperatorApplicator<AssemblyGridView,
                                                         SourceVector,
                                                         source_range_dim,
                                                         source_range_dim_cols,
                                                         SourceField,
                                                         SourceGridView,
                                                         range_range_dim,
                                                         range_range_dim_cols,
                                                         RangeField,
                                                         RangeGridView,
                                                         RangeVector>;
  using BaseType = XT::Grid::Walker<AssemblyGridView>;

public:
  using AGV = AssemblyGridView;
  using AssemblyGridViewType = AssemblyGridView;

  using SV = SourceVector;
  using SGV = SourceGridView;
  using E = XT::Grid::extract_entity_t<SGV>;
  static const constexpr size_t s_r = source_range_dim;
  static const constexpr size_t s_rC = source_range_dim_cols;
  using SF = SourceField;
  using DiscreteSourceType = ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>;
  using SourceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, SF>;

  using RV = RangeVector;
  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_range_dim;
  static const constexpr size_t r_rC = range_range_dim_cols;
  using RF = RangeField;
  using RangeType = DiscreteFunction<RV, RGV, r_r, r_rC, RF>;

  using typename BaseType::I;
  using ElementFilterType = XT::Grid::ElementFilter<AssemblyGridViewType>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AssemblyGridViewType>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AssemblyGridViewType>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AssemblyGridViewType>;

  using GenericLocalElementOperatorType = GenericLocalElementOperator<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>;
  using GenericLocalElementFunctionType = typename GenericLocalElementOperatorType::GenericFunctionType;
  using GenericLocalIntersectionOperatorType =
      GenericLocalIntersectionOperator<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>;
  using GenericLocalIntersectionFunctionType = typename GenericLocalIntersectionOperatorType::GenericFunctionType;

  LocalizableDiscreteOperatorApplicator(AssemblyGridViewType assembly_grid_view, const SourceType& src, RangeType& rng)
    : BaseType(assembly_grid_view)
    , source_(src)
    , range_(rng)
    , assembled_(false)
  {
    // to detect assembly
    this->append(
        [](/*prepare nothing*/) {}, [](const auto&) { /*apply nothing*/ }, [&](/*finalize*/) { assembled_ = true; });
  }

  const SourceType& source() const
  {
    return source_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  RangeType& range()
  {
    return range_;
  }

  using BaseType::append;

  ThisType& append(const LocalElementOperatorInterface<SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>& local_operator,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    this->append(make_local_element_operator_applicator(local_operator, range_, param).release(), filter);
    return *this;
  }

  ThisType& append(GenericLocalElementFunctionType generic_function,
                   const XT::Common::Parameter& param = {},
                   const ElementFilterType& filter = ApplyOnAllElements())
  {
    this->append(GenericLocalElementOperatorType(source(), generic_function, param.type()), param, filter);
    return *this;
  }

  ThisType&
  append(const LocalIntersectionOperatorInterface<I, SV, SGV, s_r, s_rC, SF, r_r, r_rC, RF, RGV, RV>& local_operator,
         const XT::Common::Parameter& param = {},
         const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    this->append(make_local_intersection_operator_applicator(local_operator, range_, param).release(), filter);
    return *this;
  }

  ThisType& append(GenericLocalIntersectionFunctionType generic_function,
                   const XT::Common::Parameter& param = {},
                   const IntersectionFilterType& filter = ApplyOnAllIntersections())
  {
    this->append(GenericLocalIntersectionOperatorType(source(), generic_function, param.type()), param, filter);
    return *this;
  }

  void assemble(const bool use_tbb = false)
  {
    if (assembled_)
      return;
    // This clears all appended operators, which is ok, since we are done after assembling once!
    this->walk(use_tbb);
    assembled_ = true;
  }

protected:
  const SourceType& source_;
  RangeType& range_;
  bool assembled_;
}; // class LocalizableDiscreteOperatorApplicator


template <class AssemblyGridView,
          class SourceVector,
          size_t source_range_dim = 1,
          size_t source_range_dim_cols = 1,
          class SourceField = double,
          class SourceGridView = AssemblyGridView,
          size_t range_range_dim = source_range_dim,
          size_t range_range_dim_cols = source_range_dim_cols,
          class RangeField = SourceField,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
using LocalizableOperatorBase DXT_DEPRECATED_MSG("Use LocalizableDiscreteOperatorApplicator instead (12.09.2019)!") =
    LocalizableDiscreteOperatorApplicator<AssemblyGridView,
                                          SourceVector,
                                          source_range_dim,
                                          source_range_dim_cols,
                                          SourceField,
                                          SourceGridView,
                                          range_range_dim,
                                          range_range_dim_cols,
                                          RangeField,
                                          RangeGridView,
                                          RangeVector>;


template <class AGV,
          class SV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class SGV,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV>
std::enable_if_t<XT::Grid::is_layer<AGV>::value,
                 LocalizableDiscreteOperatorApplicator<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>
make_localizable_operator_applicator(AGV assembly_grid_view,
                                     const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                                     DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range)
{
  return LocalizableDiscreteOperatorApplicator<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>(
      assembly_grid_view, source, range);
}

template <class AGV,
          size_t s_r,
          size_t s_rC,
          class SF,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV,
          class SV = RV>
std::enable_if_t<XT::Grid::is_layer<AGV>::value,
                 LocalizableDiscreteOperatorApplicator<AGV, SV, s_r, s_rC, SF, AGV, r_r, r_rC, RF, RGV, RV>>
make_localizable_operator_applicator(
    AGV assembly_grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<AGV>, s_r, s_rC, SF>& source,
    DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range)
{
  return LocalizableDiscreteOperatorApplicator<AGV, SV, s_r, s_rC, SF, AGV, r_r, r_rC, RF, RGV, RV>(
      assembly_grid_view, source, range);
}

template <class AGV,
          class SV,
          size_t s_r,
          size_t s_rC,
          class SF,
          class SGV,
          size_t r_r,
          size_t r_rC,
          class RF,
          class RGV,
          class RV>
std::enable_if_t<XT::Grid::is_layer<AGV>::value,
                 LocalizableOperatorBase<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>>
make_localizable_operator(
    AGV assembly_grid_view,
    const XT::Functions::GridFunctionInterface<XT::Grid::extract_entity_t<SGV>, s_r, s_rC, SF>& source,
    DiscreteFunction<RV, RGV, r_r, r_rC, RF>& range)
{
  return LocalizableOperatorBase<AGV, SV, s_r, s_rC, SF, SGV, r_r, r_rC, RF, RGV, RV>(
      assembly_grid_view, source, range);
}


/**
 * \note See OperatorInterface for a description of the template arguments.
 *
 * \sa OperatorInterface
 */
template <class M,
          class AssemblyGridView,
          size_t s = 1,
          size_t sC = 1,
          size_t r = s,
          size_t rC = sC,
          class RGV = AssemblyGridView,
          class SGV = AssemblyGridView>
class LocalizableOperator : public OperatorInterface<M, SGV, s, sC, r, rC, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  using ThisType = LocalizableOperator;
  using BaseType = OperatorInterface<M, SGV, s, sC, r, rC, RGV>;

public:
  using AGV = AssemblyGridView;
  using typename BaseType::F;
  using typename BaseType::V;

  using I = XT::Grid::extract_intersection_t<SGV>;

  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionInterfaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using LocalElementOperatorType = LocalElementOperatorInterface<V, SGV, s, sC, F, r, rC, F, RGV, V>;
  using LocalIntersectionOperatorType = LocalIntersectionOperatorInterface<I, V, SGV, s, sC, F, r, rC, F, RGV, V>;

  LocalizableOperator(const AGV& assembly_grid_view,
                      const SourceSpaceType& source_space,
                      const RangeSpaceType& range_space,
                      const bool linear = true,
                      const bool use_tbb = false)
    : BaseType()
    , assembly_grid_view_(assembly_grid_view)
    , source_space_(source_space)
    , range_space_(range_space)
    , linear_(linear)
    , use_tbb_(use_tbb)
  {}

  LocalizableOperator(ThisType&& source) = default;

  bool linear() const override final
  {
    return linear_;
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  ThisType& append(const LocalElementOperatorType& local_operator,
                   const XT::Grid::ElementFilter<AGV>& filter = XT::Grid::ApplyOn::AllElements<AGV>())
  {
    linear_ = linear_ && local_operator.linear();
    this->extend_parameter_type(local_operator.parameter_type());
    local_element_operators_.emplace_back(local_operator.copy(), filter.copy());
    return *this;
  }

  ThisType& append(const LocalIntersectionOperatorType& local_operator,
                   const XT::Grid::IntersectionFilter<AGV>& filter = XT::Grid::ApplyOn::AllIntersections<AGV>())
  {
    linear_ = linear_ && local_operator.linear();
    this->extend_parameter_type(local_operator.parameter_type());
    local_intersection_operators_.emplace_back(local_operator.copy(), filter.copy());
    return *this;
  }

  using BaseType::apply;

  void apply(const SourceFunctionInterfaceType& source_function,
             VectorType& range,
             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    range.set_all(0);
    auto range_function = make_discrete_function(this->range_space_, range);
    // set up the actual operator
    auto localizable_op =
        make_localizable_operator_applicator(this->assembly_grid_view_, source_function, range_function);
    // - element contributions
    for (const auto& op_and_filter : local_element_operators_) {
      const auto local_op = op_and_filter.first->with_source(source_function);
      const auto& filter = *op_and_filter.second;
      localizable_op.append(*local_op, param, filter);
    }
    // - intersection contributions
    for (const auto& op_and_filter : local_intersection_operators_) {
      const auto local_op = op_and_filter.first->with_source(source_function);
      const auto& filter = *op_and_filter.second;
      localizable_op.append(*local_op, param, filter);
    }
    // and apply it in a grid walk
    localizable_op.assemble(use_tbb_);
    DUNE_THROW_IF(!range.valid(), Exceptions::operator_error, "range contains inf or nan!");
  } // ... apply(...)

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param = {}) const override
  {
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    const auto source_function = make_discrete_function(this->source_space_, source);
    apply(source_function, range, param);
  } // ... apply(...)

  std::vector<std::string> jacobian_options() const override final
  {
    return {"finite-differences"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->jacobian_options().at(0), Exceptions::operator_error, "type = " << type);
    return {{"type", type}, {"eps", "1e-7"}};
  }

  using BaseType::jacobian;

  void jacobian(const VectorType& source,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override final
  {
    // some checks
    DUNE_THROW_IF(!source.valid(), Exceptions::operator_error, "source contains inf or nan!");
    DUNE_THROW_IF(!(this->parameter_type() <= param.type()),
                  Exceptions::operator_error,
                  "this->parameter_type() = " << this->parameter_type() << "\n   param.type() = " << param.type());
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, opts);
    DUNE_THROW_IF(opts.get<std::string>("type") != jacobian_options().at(0), Exceptions::operator_error, opts);
    const auto default_opts = jacobian_options(jacobian_options().at(0));
    const auto eps = opts.get("eps", default_opts.template get<double>("eps"));
    const auto parameter = param + XT::Common::Parameter({"finite-difference-jacobians.eps", eps});
    // append the same local ops with the same filters as in apply() above
    // - element contributions
    const auto source_function = make_discrete_function(this->source_space_, source);
    for (const auto& op_and_filter : local_element_operators_) {
      const auto local_op = op_and_filter.first->with_source(source_function);
      const auto& filter = *op_and_filter.second;
      jacobian_op.append(*local_op, source, parameter, filter);
    }
    // - intersection contributions
    for (const auto& op_and_filter : local_intersection_operators_) {
      const auto local_op = op_and_filter.first->with_source(source_function);
      const auto& filter = *op_and_filter.second;
      jacobian_op.append(*local_op, source, parameter, filter);
    }
  } // ... jacobian(...)

protected:
  const AGV assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  bool linear_;
  const bool use_tbb_;
  std::list<std::pair<std::unique_ptr<LocalElementOperatorType>, std::unique_ptr<XT::Grid::ElementFilter<AGV>>>>
      local_element_operators_;
  std::list<
      std::pair<std::unique_ptr<LocalIntersectionOperatorType>, std::unique_ptr<XT::Grid::IntersectionFilter<AGV>>>>
      local_intersection_operators_;
}; // class LocalizableOperator


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LOCALIZABLE_OPERATOR_HH
