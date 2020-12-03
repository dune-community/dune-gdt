// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2020)

#ifndef DUNE_GDT_OPERATORS_OPERATOR_HH
#define DUNE_GDT_OPERATORS_OPERATOR_HH

#include <dune/gdt/local/assembler/operator-fd-jacobian-assemblers.hh>

#include "interfaces.hh"
#include "forward-operator.hh"

namespace Dune {
namespace GDT {


template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AGV,
          class RGV = AGV>
class Operator : public OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
{
  using ThisType = Operator;
  using BaseType = OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

public:
  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::DiscreteSourceFunctionType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::V;
  using typename BaseType::VectorType;

  using I = XT::Grid::extract_intersection_t<AGV>;

  using ElementFilterType = XT::Grid::ElementFilter<AGV>;
  using IntersectionFilterType = XT::Grid::IntersectionFilter<AGV>;
  using ApplyOnAllElements = XT::Grid::ApplyOn::AllElements<AGV>;
  using ApplyOnAllIntersections = XT::Grid::ApplyOn::AllIntersections<AGV>;

  using LocalElementOperatorType = LocalElementOperatorInterface<V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV, V>;
  using LocalIntersectionOperatorType =
      LocalIntersectionOperatorInterface<I, V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV, V>;

  Operator(const AssemblyGridViewType& assembly_grid_vw,
           const SourceSpaceType& source_spc,
           const RangeSpaceType& range_spc,
           const bool requires_assembly = false,
           const std::string& logging_prefix = "",
           const std::array<bool, 3>& logging_state = {{false, false, true}})
    : BaseType({}, logging_prefix.empty() ? "Operator" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_vw)
    , source_space_(source_spc)
    , range_space_(range_spc)
    , requires_assembly_(requires_assembly)
    , linear_(true)
  {}

  Operator(const ThisType& other)
    : BaseType(other)
    , assembly_grid_view_(other.assembly_grid_view_)
    , source_space_(other.source_space_)
    , range_space_(other.range_space_)
    , requires_assembly_(other.requires_assembly_)
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
  } // Operator(...)

  Operator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;
  using BaseType::jacobian;

  /// \brief allows to fix the arguments to apply(), the resulting assembler can be appended to a GridWalker
  auto with(SourceFunctionType source_function, VectorType& range_vector, const XT::Common::Parameter& param = {}) const
  {
    this->assert_matching_range(range_vector);
    return ForwardOperatorAssembler<ThisType>(
        *this, source_function, range_vector, param, this->logger.prefix + "_assembler", this->logger.state);
  }

  /// \brief allows to fix the arguments to apply(), the resulting assembler can be appended to a GridWalker
  auto with(VectorType& source_vector, VectorType& range_vector, const XT::Common::Parameter& param = {}) const
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    return ForwardOperatorAssembler<ThisType>(*this,
                                              DiscreteSourceFunctionType(source_space_, source_vector),
                                              range_vector,
                                              param,
                                              this->logger.prefix + "_assembler",
                                              this->logger.state);
  } // ... with(...)

  /// \name Required by OperatorInterface
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
    DUNE_THROW_IF(requires_assembly_,
                  Exceptions::operator_error,
                  "You need to call assemble() first (or a derived class forgot to set requires_assembly_ to true in "
                  "assemble())!");
    this->assert_matching_range(range_vector);
    auto assembler = this->with(source_function, range_vector, param);
    XT::Grid::Walker<AGV> walker(this->assembly_grid_view_);
    walker.append(assembler);
    walker.walk(/*use_tbb=*/true);
  } // ... apply(...)

  /// \}
  /// \name Required by OperatorInterface
  /// \{

  const AssemblyGridViewType& assembly_grid_view() const override
  {
    return assembly_grid_view_;
  }

  const SourceSpaceType& source_space() const override
  {
    return source_space_;
  }

  /// \brief Simply redirects to apply(source_function, range_vector) to avoid the base interpolation + redirect
  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override
  {
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    this->apply(make_discrete_function(source_space_, source_vector), range_vector, param);
  }

protected:
  std::vector<XT::Common::Configuration> all_jacobian_options() const override
  {
    return {{{"type", "finite-differences"}, {"eps", "1e-7"}}};
  }

public:
  void jacobian(const VectorType& source_vector,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.sup_norm()=" << jacobian_op.matrix().sup_norm() << print(opts, {{"oneline", "true"}})
                << ", param=" << param << ")" << std::endl;
    this->assert_jacobian_opts(opts); // ensures that type finite-differences is requested
    this->assert_matching_source(source_vector);
    const std::string type = opts.get<std::string>("type");
    const auto default_opts = this->jacobian_options(type);
    const double eps = opts.get("eps", default_opts.template get<double>("eps"));
    const auto parameter = param + XT::Common::Parameter({"finite-difference-jacobians.eps", eps});
    // append the same local ops with the same filters as in apply() above
    LOG_(info) << "appending {" << element_data_.size() << "|" << intersection_data_.size()
               << "} local {element|intersection} operators to jacobian_op ..." << std::endl;
    const auto source_function = make_discrete_function(source_space_, source_vector);
    const auto append_local_operator_to_jacobian = [&](const auto& list_of_data) {
      for (const auto& data : list_of_data) {
        const auto local_operator = data.first->with_source(source_function);
        const auto& filter = *data.second;
        jacobian_op.append(*local_operator, source_vector, parameter, filter);
      }
    };
    append_local_operator_to_jacobian(element_data_);
    append_local_operator_to_jacobian(intersection_data_);
  } // ... jacobian(...)

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
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  bool requires_assembly_;
  bool linear_;
  std::list<std::pair<std::unique_ptr<LocalElementOperatorType>, std::unique_ptr<ElementFilterType>>> element_data_;
  std::list<std::pair<std::unique_ptr<LocalIntersectionOperatorType>, std::unique_ptr<IntersectionFilterType>>>
      intersection_data_;
}; // class Operator


template <class MatrixType, // <- needs to be manually specified
          class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC>
auto make_operator(const AssemblyGridViewType& assembly_grid_view,
                   const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                   const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  return Operator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}

template <class AssemblyGridViewType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
auto make_operator(const AssemblyGridViewType& assembly_grid_view,
                   const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                   const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  return make_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}


template <class MatrixType, // <- needs to be manually specified
          class GV,
          size_t r,
          size_t rC,
          class F>
auto make_operator(const SpaceInterface<GV, r, rC, F>& space,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return Operator<GV, r, rC, r, rC, F, MatrixType, GV, GV>(
      space.grid_view(), space, space, logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F>
auto make_operator(const SpaceInterface<GV, r, rC, F>& space,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  return make_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(space, logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_OPERATOR_HH
