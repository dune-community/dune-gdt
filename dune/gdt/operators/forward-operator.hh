// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_OPERATORS_FORWARD_OPERATOR_HH
#define DUNE_GDT_OPERATORS_FORWARD_OPERATOR_HH

#include "interfaces.hh"

namespace Dune {
namespace GDT {


// forward, required in ForwardOperator
template <class ForwardOperatorType>
class ForwardOperatorAssembler;


template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class Operator; // include is below


template <class AssemblyGridView,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class V = XT::LA::IstlDenseVector<F>,
          class SGV = AssemblyGridView,
          class RGV = AssemblyGridView>
class ForwardOperator : public ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, V, RGV>
{
  using ThisType = ForwardOperator;
  using BaseType = ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, V, RGV>;

public:
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::VectorType;

  using AssemblyGridViewType = AssemblyGridView;
  using AGV = AssemblyGridView;
  using I = XT::Grid::extract_intersection_t<AGV>;

  using ElementFilterType = XT::Grid::ElementFilter<AGV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AGV>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AGV>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AGV>;

  using LocalElementOperatorType = LocalElementOperatorInterface<V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV, V>;
  using LocalIntersectionOperatorType =
      LocalIntersectionOperatorInterface<I, V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV, V>;

  ForwardOperator(const AssemblyGridViewType& assembly_grid_vw,
                  const RangeSpaceType& range_spc,
                  const std::string& logging_prefix = "",
                  const std::array<bool, 3>& logging_state = {{false, false, true}})
    : BaseType({}, logging_prefix.empty() ? "ForwardOperator" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_vw)
    , range_space_(range_spc)
    , linear_(true)
  {}

  ForwardOperator(const ThisType& other)
    : BaseType(other)
    , assembly_grid_view_(other.assembly_grid_view_)
    , range_space_(other.range_space_)
    , linear_(other.linear_)
  {
    const auto copy_local_data = [](const auto& origin, auto& target) {
      for (const auto& data : origin) {
        const auto& local_operator = *data.first;
        const auto& filter = *data.second;
        target.emplace_back(local_operator.copy(), filter.copy());
      }
    };
    copy_local_data(other.element_data_, element_data_);
    copy_local_data(other.intersection_data_, intersection_data_);
  } // ForwardOperator(...)

  ForwardOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;

  const AssemblyGridViewType& assembly_grid_view() const
  {
    return assembly_grid_view_;
  }

  /// \brief allows to fix the arguments to apply(), the resulting assembler can be appended to a GridWalker
  auto with(SourceFunctionType source_function, VectorType& range_vector, const XT::Common::Parameter& param = {}) const
  {
    this->assert_matching_range(range_vector);
    return ForwardOperatorAssembler<ThisType>(
        *this, source_function, range_vector, param, this->logger.prefix + "_assembler", this->logger.state);
  }

  /// \brief allows to obtain the corresponding Operator
  template <class MatrixType = XT::LA::matrix_t<V>>
  auto with(const SpaceInterface<SGV, s_r, s_rC, F>& source_space) const
  {
    static_assert(XT::LA::is_matrix<MatrixType>::value, "");
    Operator<AGV, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV> discrete_op(
        assembly_grid_view_, source_space, range_space_, this->logger.prefix, this->logger.state);
    const auto append_data = [&](const auto& origin) {
      for (const auto& data : origin) {
        const auto& local_operator = *data.first;
        const auto& filter = *data.second;
        discrete_op += {local_operator, filter};
      }
    };
    append_data(element_data_);
    append_data(intersection_data_);
    return discrete_op;
  } // ... with(...)

  /// \name Required by ForwardOperatorInterface
  /// \{

  const RangeSpaceType& range_space() const override
  {
    return range_space_;
  }

  bool linear() const override
  {
    return linear_;
  }

  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override
  {
    this->assert_matching_range(range_vector);
    auto assembler = this->with(source_function, range_vector, param);
    XT::Grid::Walker<AGV> walker(this->assembly_grid_view_);
    walker.append(assembler);
    walker.walk(/*use_tbb=*/true);
  }

  /// \}
  /// \name These methods allow to append local operators
  ///\{

  ThisType& operator+=(const LocalElementOperatorType& local_operator)
  {
    LOG_(info) << "+=(local_element_operator=" << &local_operator << ")" << std::endl;
    this->extend_parameter_type(local_operator.parameter_type());
    linear_ = linear_ && local_operator.linear();
    element_data_.emplace_back(local_operator.copy(), new ApplyOnAllElements());
    return *this;
  }

  ThisType&
  operator+=(const std::pair<const LocalElementOperatorType&, const ElementFilterType&>& local_operator__filter)
  {
    const auto& local_operator = local_operator__filter.first;
    const auto& filter = local_operator__filter.second;
    LOG_(info) << "+=(local_element_operator=" << &local_operator << ", element_filter=" << &filter << ")" << std::endl;
    this->extend_parameter_type(local_operator.parameter_type());
    linear_ = linear_ && local_operator.linear();
    element_data_.emplace_back(local_operator.copy(), filter.copy());
    return *this;
  } // ... operator+=(...)

  ThisType& operator+=(const LocalIntersectionOperatorType& local_operator)
  {
    LOG_(info) << "+=(local_intersection_operator=" << &local_operator << ")" << std::endl;
    this->extend_parameter_type(local_operator.parameter_type());
    linear_ = linear_ && local_operator.linear();
    intersection_data_.emplace_back(local_operator.copy(), new ApplyOnAllIntersections());
    return *this;
  }

  ThisType& operator+=(
      const std::pair<const LocalIntersectionOperatorType&, const IntersectionFilterType&>& local_operator__filter)
  {
    const auto& local_operator = local_operator__filter.first;
    const auto& filter = local_operator__filter.second;
    LOG_(info) << "+=(local_element_operator=" << &local_operator << ", intersection_filter=" << &filter << ")"
               << std::endl;
    this->extend_parameter_type(local_operator.parameter_type());
    linear_ = linear_ && local_operator.linear();
    intersection_data_.emplace_back(local_operator.copy(), filter.copy());
    return *this;
  } // ... operator+=(...)

  /// \}

  const std::list<std::pair<std::unique_ptr<LocalElementOperatorType>, std::unique_ptr<ElementFilterType>>>&
  element_data() const
  {
    return element_data_;
  }

  const std::list<std::pair<std::unique_ptr<LocalIntersectionOperatorType>, std::unique_ptr<IntersectionFilterType>>>&
  intersection_data() const
  {
    return intersection_data_;
  }

protected:
  const AssemblyGridViewType& assembly_grid_view_;
  const RangeSpaceType& range_space_;
  bool linear_;
  std::list<std::pair<std::unique_ptr<LocalElementOperatorType>, std::unique_ptr<ElementFilterType>>> element_data_;
  std::list<std::pair<std::unique_ptr<LocalIntersectionOperatorType>, std::unique_ptr<IntersectionFilterType>>>
      intersection_data_;
}; // class ForwardOperator


/// \note This assembler can be used for ForwardOperator and Operator, we thus use duck-typing
template <class ForwardOperatorType>
class ForwardOperatorAssembler : public XT::Grid::ElementAndIntersectionFunctor<typename ForwardOperatorType::AGV>
{
  using ThisType = ForwardOperatorAssembler;
  using BaseType = XT::Grid::ElementAndIntersectionFunctor<typename ForwardOperatorType::AGV>;

public:
  using typename BaseType::E;
  using typename BaseType::I;

  static constexpr size_t s_r = ForwardOperatorType::s_r;
  static constexpr size_t s_rC = ForwardOperatorType::s_rC;
  static constexpr size_t r_r = ForwardOperatorType::r_r;
  static constexpr size_t r_rC = ForwardOperatorType::r_rC;
  using F = typename ForwardOperatorType::F;
  using DiscreteRangeFunctionType = typename ForwardOperatorType::DiscreteRangeFunctionType;
  using SourceFunctionType = typename ForwardOperatorType::SourceFunctionType;
  using VectorType = typename ForwardOperatorType::VectorType;
  using LocalElementOperatorType = typename ForwardOperatorType::LocalElementOperatorType;
  using LocalIntersectionOperatorType = typename ForwardOperatorType::LocalIntersectionOperatorType;
  using ElementFilterType = typename ForwardOperatorType::ElementFilterType;
  using IntersectionFilterType = typename ForwardOperatorType::IntersectionFilterType;

  ForwardOperatorAssembler(ForwardOperatorType oprtr,
                           SourceFunctionType src,
                           VectorType& range_vector,
                           const XT::Common::Parameter& param = {},
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = {{false, false, true}})
    : BaseType(logging_prefix.empty() ? "ForwardOperatorAssembler" : logging_prefix, logging_state)
    , operator_(oprtr)
    , source_(src.copy_as_grid_function())
    , range_vector_(range_vector)
    , param_(param)
    , range_function_(operator_.range_space(), range_vector_)
    , local_range_inside_(range_function_.local_discrete_function())
    , local_range_outside_(range_function_.local_discrete_function())
  {}

  ForwardOperatorAssembler(const ThisType& other)
    : operator_(other.operator_)
    , source_(other.source_->copy_as_grid_function())
    , range_vector_(other.range_vector_)
    , param_(other.param_)
    , range_function_(operator_.range_space(), range_vector_)
    , local_range_inside_(range_function_.local_discrete_function())
    , local_range_outside_(range_function_.local_discrete_function())
  {}

  ForwardOperatorAssembler(ThisType&& source) = default;

  /// \name Required by ElementAndIntersectionFunctor.
  /// \{

  BaseType* copy() override
  {
    return new ThisType(*this);
  }

  void prepare() override
  {
    // clear range
    range_vector_.set_all(0);
    // set source in local operators
    const auto set_source_in_local_ops = [&](const auto& origin, auto& target) {
      for (const auto& data : origin) {
        const auto& local_operator = *data.first;
        const auto& filter = *data.second;
        target.emplace_back(local_operator.with_source(*source_), filter.copy());
      }
    };
    set_source_in_local_ops(operator_.element_data(), element_data_);
    set_source_in_local_ops(operator_.intersection_data(), intersection_data_);
  } // ... prepare(...)

  void apply_local(const E& element) override
  {
    if (element_data_.size() > 0)
      local_range_inside_->bind(element);
    for (auto& data : element_data_) {
      auto& local_operator = *data.first;
      const auto& filter = *data.second;
      if (filter.contains(operator_.assembly_grid_view(), element)) {
        local_operator.bind(element);
        local_operator.apply(*local_range_inside_, param_);
      }
    }
  } // ... apply_local(...)

  void apply_local(const I& intersection, const E& inside_element, const E& outside_element) override
  {
    if (intersection_data_.size() > 0) {
      local_range_inside_->bind(inside_element);
      local_range_outside_->bind(outside_element);
    }
    for (auto& data : intersection_data_) {
      auto& local_operator = *data.first;
      const auto& filter = *data.second;
      if (filter.contains(operator_.assembly_grid_view(), intersection)) {
        local_operator.bind(intersection);
        local_operator.apply(*local_range_inside_, *local_range_outside_, param_);
      }
    }
  } // ... apply_local(...

  void finalize() override
  {
    element_data_.clear();
    intersection_data_.clear();
  }

  /// \}

protected:
  const ForwardOperatorType operator_;
  const std::unique_ptr<XT::Functions::GridFunctionInterface<E, s_r, s_rC, F>> source_;
  VectorType& range_vector_;
  const XT::Common::Parameter param_;
  DiscreteRangeFunctionType range_function_;
  std::unique_ptr<typename DiscreteRangeFunctionType::LocalDiscreteFunctionType> local_range_inside_;
  std::unique_ptr<typename DiscreteRangeFunctionType::LocalDiscreteFunctionType> local_range_outside_;
  std::list<std::pair<std::unique_ptr<LocalElementOperatorType>, std::unique_ptr<ElementFilterType>>> element_data_;
  std::list<std::pair<std::unique_ptr<LocalIntersectionOperatorType>, std::unique_ptr<IntersectionFilterType>>>
      intersection_data_;
}; // class ForwardOperatorAssembler


template <class AssemblyGridViewType,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class F,
          size_t s_r = r_r,
          size_t s_rC = r_rC,
          class V = XT::LA::IstlDenseVector<F>,
          class SGV = RGV>
auto make_forward_operator(const AssemblyGridViewType& assembly_grid_view,
                           const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  return ForwardOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, V, SGV, RGV>(
      assembly_grid_view, range_space, logging_prefix, logging_state);
}


template <class GV,
          size_t r_r,
          size_t r_rC,
          class F,
          size_t s_r = r_r,
          size_t s_rC = r_rC,
          class V = XT::LA::IstlDenseVector<F>,
          class SGV = GV>
auto make_forward_operator(const SpaceInterface<GV, r_r, r_rC, F>& space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  return ForwardOperator<GV, s_r, s_rC, r_r, r_rC, F, V, SGV, GV>(
      space.grid_view(), space, logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#include "operator.hh"

#endif // DUNE_GDT_OPERATORS_FORWARD_OPERATOR_HH
