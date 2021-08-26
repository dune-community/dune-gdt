// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2018)
//   René Fritze     (2016, 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_INTERFACES_HH
#define DUNE_GDT_OPERATORS_INTERFACES_HH

#include <cmath>
#include <type_traits>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/functions/grid-function.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/solver.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/algorithms/newton.hh>
#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/print.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

namespace Dune {
namespace GDT {


// forwards
// - required for numeric operators, includes are below
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class ConstantOperator;

template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class ConstLincombOperator;

template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class LincombOperator;

// - required for the jacobian, includes are below
template <class AGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class F, class M, class SGV, class RGV>
class MatrixOperator;


template <class SourceGridView,
          size_t source_dim = 1,
          size_t source_dim_cols = 1,
          size_t range_dim = source_dim,
          size_t range_dim_cols = source_dim_cols,
          class Field = double,
          class RangeGridView = SourceGridView>
class BilinearFormInterface
  : public XT::Common::ParametricInterface
  , public XT::Common::WithLogger<BilinearFormInterface<SourceGridView,
                                                        source_dim,
                                                        source_dim_cols,
                                                        range_dim,
                                                        range_dim_cols,
                                                        Field,
                                                        RangeGridView>>
{
  static_assert(XT::Grid::is_view<SourceGridView>::value);
  static_assert(XT::Grid::is_view<RangeGridView>::value);
  static_assert(SourceGridView::dimension == RangeGridView::dimension);
  using ThisType = BilinearFormInterface;
  using Logger = XT::Common::WithLogger<BilinearFormInterface<SourceGridView,
                                                              source_dim,
                                                              source_dim_cols,
                                                              range_dim,
                                                              range_dim_cols,
                                                              Field,
                                                              RangeGridView>>;

public:
  using FieldType = Field;
  using F = FieldType;

  using SGV = SourceGridView;
  using SE = XT::Grid::extract_entity_t<SGV>;
  static constexpr size_t s_r = source_dim;
  static constexpr size_t s_rC = source_dim_cols;

  using SourceFunctionType = XT::Functions::GridFunction<SE, s_r, s_rC, F>;

  using RGV = RangeGridView;
  using RE = XT::Grid::extract_entity_t<SGV>;
  static constexpr size_t r_r = range_dim;
  static constexpr size_t r_rC = range_dim_cols;

  using RangeFunctionType = XT::Functions::GridFunction<RE, r_r, r_rC, F>;

  explicit BilinearFormInterface(const XT::Common::ParameterType& param_type = {},
                                 const std::string& logging_prefix = "",
                                 const std::array<bool, 3>& logging_enabled = XT::Common::default_logger_state())
    : XT::Common::ParametricInterface(param_type)
    , Logger(logging_prefix.empty() ? "BilinearFormInterface" : logging_prefix, logging_enabled)
  {
    LOG_(debug) << "BilinearFormInterface(param_type=" << print(param_type) << ")" << std::endl;
  }

  BilinearFormInterface(const ThisType& other) = default;

  BilinearFormInterface(ThisType&& source) = default;

  virtual ~BilinearFormInterface() = default;

  /// \name These methods have to be implemented to define the action of the BilinearForm.
  /// \{

  virtual FieldType apply2(RangeFunctionType /*range_function*/,
                           SourceFunctionType /*source_function*/,
                           const XT::Common::Parameter& /*param*/ = {}) const = 0;

  /// \}

  /// \brief Allows the implementation to do preparatory work (i.e., assemble the matrix of a matrix-based operator).
  /// \note In general, you have to call this method before calling apply2!
  /// \todo Drop use_tbb, has to be specified on construction or somewhere else
  virtual void assemble(const bool use_tbb = false)
  {
    LOG_(debug) << "assemble(use_tbb=" << use_tbb << ")" << std::endl;
    LOG_(info) << "not assembling bilinear form" << std::endl;
  }

  FieldType norm(RangeFunctionType range_function, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "norm(range_function=" << &range_function << ", param=" << print(param) << ")" << std::endl;
    if constexpr (std::is_same_v<RE, SE> && (r_r == s_r) && (r_rC == s_rC)) {
      LOG_(debug) << "calling apply2(range_function, range_function, param)... " << std::endl;
      const auto result = this->apply2(range_function, range_function, param);
      if (XT::Common::is_zero(result)) {
        LOG_(info) << "result of apply2() is nearly zero (" << result << "), returning 0" << std::endl;
        return FieldType(0);
      }
      DUNE_THROW_IF(result < 0,
                    Exceptions::operator_error,
                    "apply2(range, range) of BilinearForm did not yield a positive value (" << result << ")!");
      return std::sqrt(result);
    } else {
      DUNE_THROW(Exceptions::operator_error, "BilinearForm does not induce a norm, source and range differ!");
      return FieldType(0);
    }
  } // ... norm(...)
}; // class BilinearFormInterface


template <class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class Vector = XT::LA::IstlDenseVector<F>,
          class RGV = SGV>
class ForwardOperatorInterface : public BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

public:
  using ThisType = ForwardOperatorInterface;
  using BaseType = BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>;

  using typename BaseType::FieldType;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::SourceFunctionType;

  using VectorType = Vector;
  using V = VectorType;

  using RangeSpaceType = SpaceInterface<RGV, r_r, r_rC, F>;
  using DiscreteRangeFunctionType = DiscreteFunction<V, RGV, r_r, r_rC, F>;

  explicit ForwardOperatorInterface(const XT::Common::ParameterType& param_type = {},
                                    const std::string& logging_prefix = "",
                                    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType(param_type, logging_prefix.empty() ? "ForwardOperatorInterface" : logging_prefix, logging_state)
  {
    LOG_(debug) << "ForwardOperatorInterface(param_type=" << print(param_type) << ")" << std::endl;
  }

  ForwardOperatorInterface(const ThisType& other) = default;

  ForwardOperatorInterface(ThisType&& source) = default;

  virtual ~ForwardOperatorInterface() = default;

  // pull in methods from BilinearFormInterface
  using BaseType::apply2;

  /// \name These methods have to be implemented to define the discrete range of the Operator.
  /// \{

  virtual const RangeSpaceType& range_space() const = 0;

  /// \}
  /// \name These methods have to be implemented to define the forward-action of the Operator.
  /// \{

  virtual bool linear() const = 0;

  virtual void apply(SourceFunctionType source_function,
                     VectorType& range_vector,
                     const XT::Common::Parameter& param = {}) const = 0;

  /// \}
  /// \name These methods are required by BilinearFormInterface and default-implemented here.
  /// \{

  /// \brief Carries out an interpolation of range_function and calls the other apply2() afterwards (not optimal!).
  virtual FieldType apply2(RangeFunctionType range_function,
                           SourceFunctionType source_function,
                           const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "apply2(range_function=" << &range_function << ", source_function=" << &source_function
                << ", param=" << print(param) << ")" << std::endl;
    LOG_(info) << "interpolating range_function ..." << std::endl;
    auto discrete_range_function = interpolate<V>(range_function, this->range_space(), param);
    LOG_(debug) << "calling discrete apply2() variant ..." << std::endl;
    return this->apply2(discrete_range_function.dofs().vector(), source_function, param);
  }

  /// \}
  /// \name These methods extend the ones from BilinearFormInterface.
  /// \{

  virtual FieldType apply2(const VectorType& range_vector,
                           SourceFunctionType source_function,
                           const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply2(range_vector.sup_norm()=" << range_vector.sup_norm()
                << ", source_function=" << &source_function << ", param=" << print(param) << ")" << std::endl;
    this->assert_matching_range(range_vector);
    LOG_(info) << "computing apply(source_function, param).dot(range_vector) ..." << std::endl;
    return this->apply(source_function, param).dofs().vector().dot(range_vector);
  }

  /// \}
  /// \name These apply methods are provided for convenience.
  /// \{

  virtual void apply(SourceFunctionType source_function,
                     DiscreteRangeFunctionType& discrete_range_function,
                     const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", discrete_range_function=" << &discrete_range_function << ", param=" << print(param) << ")"
                << "\n"
                << "  redicrecting to discrete apply() variant ..." << std::endl;
    this->apply(source_function, discrete_range_function.dofs().vector(), param);
  }

  virtual DiscreteRangeFunctionType apply(SourceFunctionType source_function,
                                          const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply(source_function=" << &source_function << ", param=" << print(param) << ")"
                << "\n"
                << "  creating discrete_range_function and redicrecting to discrete apply() variant ..." << std::endl;
    DiscreteRangeFunctionType discrete_range_function(this->range_space());
    this->apply(source_function, discrete_range_function.dofs().vector(), param);
    return discrete_range_function;
  }

  /// \}

protected:
  void assert_matching_range(const VectorType& range_vector) const
  {
    DUNE_THROW_IF(!this->range_space().contains(range_vector),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "range_space().mapper().size() = " << this->range_space().mapper().size() << "\n"
                                                     << "   range_vector.size() = " << range_vector.size());
  } // ... assert_matching_range(...)

  void assert_matching_range(const DiscreteRangeFunctionType& discrete_range_function) const
  {
    DUNE_THROW_IF(!this->range_space().contains(discrete_range_function),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "range_space().mapper().size() = " << this->range_space().mapper().size() << "\n"
                                                     << "   discrete_range_function.dofs().vector().size() = "
                                                     << discrete_range_function.dofs().vector().size() << "\n"
                                                     << "   range_space.type() = " << this->range_space().type() << "\n"
                                                     << "   discrete_range_function.space().type() = "
                                                     << discrete_range_function.space().type());
  } // ... assert_matching_range(...)
}; // class ForwardOperatorInterface


/**
 * \brief Interface for operators (and two-forms).
 *
 * Considering (discrete) spaces V_h and W_h, and a field F, this interface models
 *
 * - operators A: V_h -> W_h and
 *
 * - two-forms B: W_h x V_h -> F.
 *
 * The source of the operator V_h is the (discrete) space
 *
 *   V_h := \{ v_h: \tau_h^s -> F^{s_r \times s_rC} \}
 *
 * (modelled by SpaceInterface) of functions mapping from a partition \tau_h^s (modelled by SourceGridView) of a
 * physical domain to a (possibly vector- or matrix-valued) vector space F^{s_r \times s_rC} (modelled by Field,
 * source_dim and source_dim_cols). The range of the operator W_h (identified with its dual, since we are in the
 * discrete setting) is the (discrete) space
 *
 *   W_h := \{ w_h: \tau_h^r -> F^{r_r \times r_rC} \}
 *
 * (modelled by SpaceInterface) of functions mapping from a partition \tau_h^r (modelled by RangeGridView) of a
 * physical domain to a (possibly vector- or matrix-valued) vector space F^{r_r \times r_rC} (modelled by Field,
 * range_dim and range_dim_cols).
 *
 * The functions v_h \in V_h (and w_h \in W_h), to which the operator can be applied to, are represented by their DoF
 * vectors, which are then interpreted as discrete functions in the respective source_space or range_space.
 *
 * The appropriate vector type (derived from XT::LA::VectorInterface) is automatically deduced from the given
 * matrix type (derived from XT::LA::MatrixInterface, modelled by Matrix), as well as the underlying field.
 *
 * \note In general, one would like to have differente fields for the source vector, the range vector, the matrix and
 *       the result of apply2(). However, this is postponed in favor of fewer template arguments, until we require it.
 */
template <class AssemblyGridView,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class Matrix = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AssemblyGridView,
          class RGV = AssemblyGridView>
class OperatorInterface : public ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, XT::LA::vector_t<Matrix>, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_matrix<Matrix>::value, "");

public:
  using ThisType = OperatorInterface;
  using BaseType = ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, XT::LA::vector_t<Matrix>, RGV>;

  using typename BaseType::FieldType;
  using typename BaseType::V;
  using typename BaseType::VectorType;

  using typename BaseType::SourceFunctionType;
  using SourceSpaceType = SpaceInterface<SGV, s_r, s_rC, F>;
  using DiscreteSourceFunctionType = DiscreteFunction<V, SGV, s_r, s_rC, F>;
  using ConstDiscreteSourceFunctionType = ConstDiscreteFunction<V, SGV, s_r, s_rC, F>;

  using typename BaseType::DiscreteRangeFunctionType;
  using ConstDiscreteRangeFunctionType = ConstDiscreteFunction<V, RGV, r_r, r_rC, F>;

  using AssemblyGridViewType = AssemblyGridView;
  using AGV = AssemblyGridView;
  using MatrixType = Matrix;
  using M = Matrix;
  using MatrixOperatorType = MatrixOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

  using ConstLincombOperatorType = ConstLincombOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;
  using LincombOperatorType = LincombOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

  explicit OperatorInterface(const XT::Common::ParameterType& param_type = {},
                             const std::string& logging_prefix = "",
                             const std::array<bool, 3>& logging_enabled = XT::Common::default_logger_state())
    : BaseType(param_type, logging_prefix.empty() ? "OperatorInterface" : logging_prefix, logging_enabled)
  {
    LOG_(debug) << "OperatorInterface(param_type=" << print(param_type) << ")" << std::endl;
  }

  OperatorInterface(const ThisType& other) = default;

  OperatorInterface(ThisType&& source) = default;

  ~OperatorInterface() override = default;

  // pull in methods from various base classes
  using BaseType::apply;
  using BaseType::apply2;
  using BaseType::norm;

  /// \name These methods have to be implemented to define the discrete source of the DiscreteOperator.
  /// \{

  virtual const SourceSpaceType& source_space() const = 0;

  /// \}
  /// \name These methods are required to define the assembly region of the operator and its jacobian.
  /// \{

  virtual const AssemblyGridViewType& assembly_grid_view() const = 0;

  /// \}
  /// \name These methods are required by ForwardOperatorInterface and default-implemented here.
  /// \{

  /// \brief Carries out an interpolation of source_function and calls the other apply() afterwards (not optimal!).
  virtual void apply(SourceFunctionType source_function,
                     VectorType& range_vector,
                     const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    LOG_(info) << "interpolating source_function ..." << std::endl;
    auto discrete_source_function = interpolate<V>(source_function, this->source_space());
    LOG_(debug) << "calling discrete apply() variant ..." << std::endl;
    this->apply(discrete_source_function.dofs().vector(), range_vector, param);
  }

  /// \name These methods extend the ones from BilinearFormInterface and OperatorInterface.
  /// \{

  virtual FieldType
  apply2(const VectorType& range_vector, const VectorType& source_vector, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply2(range_vector.sup_norm()=" << range_vector.sup_norm()
                << ", source_vector.sup_norm()=" << source_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    this->assert_matching_source(source_vector);
    LOG_(info) << "computing apply(source, param).dot(range_vector) with euklidean product ..." << std::endl;
    try {
      LOG_(debug) << "therefore trying discrete apply() variant ..." << std::endl;
      return this->apply(source_vector, param).dot(range_vector);
    } catch (Exceptions::operator_error&) {
      LOG_(debug) << "failed, trying other apply() variant ..." << std::endl;
      return this->apply(make_discrete_function(this->source_space(), source_vector), param)
          .dofs()
          .vector()
          .dot(range_vector);
    }
  } // ... apply2(...)

  FieldType norm(const VectorType& range_vector, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "norm(range_vector.sup_norm=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    if constexpr ((r_r == s_r) && (r_rC == s_rC)) {
      this->assert_matching_range(range_vector);
      this->assert_matching_source(range_vector);
      LOG_(debug) << "calling norm(discrete_range_function, param)... " << std::endl;
      // sanity checks for zero are carried out in the BilinearFormInterface
      return this->norm(DiscreteRangeFunctionType(this->range_space(), range_vector), param);
    } else {
      DUNE_THROW(Exceptions::operator_error, "DiscreteOperator does not induce a norm, source and range differ!");
      return FieldType(0);
    }
  } // ... norm(...)

  /// \}
  /// \name These methods should be implemented to define the forward-action of the DiscreteOperator.
  /// \{

  virtual void apply(const VectorType& /*source_vector*/,
                     VectorType& /*range_vector*/,
                     const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error,
               "This DiscreteOperator does not support apply(source_vector, range_vector)!");
  }

  /// \}
  /// \name These methods should be implemented to define the jacobian of the DiscreteOperator.
  /// \{

  /**
   * \brief Returns a list of string identifiers for different means to compute the jacobian (in decending order based
   *        on what-we-think-best-ness), which can be passed to jacobian().
   * \sa    all_jacobian_options
   */
  virtual std::vector<std::string> jacobian_options() const
  {
    const auto all_opts = this->all_jacobian_options();
    DUNE_THROW_IF(all_opts.empty(), Exceptions::operator_error, "This DiscreteOperator does not provide a jacobian!");
    std::vector<std::string> types;
    for (const auto& opts : all_opts) {
      DUNE_THROW_IF(!opts.has_key("type"),
                    Exceptions::operator_error,
                    "This DiscreteOperator does not correctly report its jacobian options!");
      types.emplace_back(opts.template get<std::string>("type"));
    }
    return types;
  } // ... jacobian_options(...)

  /**
   * \brief For each identifier returned by jacobian_options(), returns an options dict with at least the key "type" set
   *        to type
\code
jacobian_options(some_type).get<std::string>("type") == some_type
\endcode
   *        and possibly other key/value pairs, which can be passed to jacobian().
   * \sa    all_jacobian_options
   */
  virtual XT::Common::Configuration jacobian_options(const std::string& type) const
  {
    const auto all_opts = this->all_jacobian_options();
    DUNE_THROW_IF(all_opts.empty(), Exceptions::operator_error, "This DiscreteOperator does not provide a jacobian!");
    for (const auto& opts : all_opts) {
      DUNE_THROW_IF(!opts.has_key("type"),
                    Exceptions::operator_error,
                    "This DiscreteOperator does not correctly report its jacobian options!");
      if (opts.template get<std::string>("type") == type)
        return opts;
    }
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
               "Given type is not one of jacobian_options()!"
                   << "\n"
                   << "   type = " << type << "\n"
                   << "   jacobian_options() = " << print(this->jacobian_options()));
    return XT::Common::Configuration();
  } // ... jacobian_options(...)

protected:
  /**
   * Used by derived DiscreteOperators to feed the implementations of jacobian_options().
   * Simply return a vector of configurations with each at least a "type" key.
   */
  virtual std::vector<XT::Common::Configuration> all_jacobian_options() const
  {
    return {};
  }

public:
  /**
   * \brief Either appends suitable functors to the jacobian_op (such that the jacobian of this operator is assembled
   *        additively into jacobian_op) or adds the jacobian of this operator to jacobian_op.matrix().
   *
   * \note You need to call jacobian_op.assemble() to be sure to have the jacobian fully assembled.
   **/
  virtual void jacobian(const VectorType& source_vector,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Configuration& opts,
                        const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm()
                << ",\n   opts=" << print(opts, {{"oneline", "true"}}) << ",\n   param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_source(source_vector);
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, "opts = " << opts);
    const auto type = opts.get<std::string>("type");
    if (type == "zero") {
      LOG_(info) << "not adding zero jacobian ..." << std::endl;
      return;
    } else {
      DUNE_THROW(Exceptions::operator_error,
                 "This DiscreteOperator reports to support jacobian(source_vector, opts) with any of the opts below,"
                     << "\n"
                     << "but its implementation does not override the respective method!"
                     << "\n\n"
                     << print(opts));
    }
    DUNE_THROW(Exceptions::operator_error, "This DiscreteOperator does not support jacobian(source_vector)!");
  } // ... jacobian(...)

  /// \}

  /// \name These methods should be implemented to define the inverse action of the DiscreteOperator.
  /// \{

  /**
   * \brief Returns a list of string identifiers for different means to apply the inverse (in decending order based on
   *        what-we-think-best-ness), which can be passed to apply_inverse().
   * \sa    all_invert_options
   */
  virtual std::vector<std::string> invert_options() const
  {
    const auto all_opts = this->all_invert_options();
    DUNE_THROW_IF(all_opts.empty(), Exceptions::operator_error, "This DiscreteOperator is not invertible!");
    std::vector<std::string> types;
    for (const auto& opts : all_opts) {
      DUNE_THROW_IF(!opts.has_key("type"),
                    Exceptions::operator_error,
                    "This DiscreteOperator does not correctly report its inversion options!");
      types.emplace_back(opts.template get<std::string>("type"));
    }
    return types;
  } // ... invert_options(...)

  /**
   * \brief For each identifier returned by invert_options(), returns an options dict with at least the key "type" set
   *        to type
\code
invert_options(some_type).get<std::string>("type") == some_type
\endcode
   *        and possibly other key/value pairs, which can be passed to jacobian().
   * \sa    all_invert_options
   * \todo  Allow to pass jacobian options as subcfg in newton.
   */
  virtual XT::Common::Configuration invert_options(const std::string& type) const
  {
    const auto all_opts = this->all_invert_options();
    DUNE_THROW_IF(all_opts.empty(), Exceptions::operator_error, "This DiscreteOperator is not invertible!");
    for (const auto& opts : all_opts) {
      DUNE_THROW_IF(!opts.has_key("type"),
                    Exceptions::operator_error,
                    "This DiscreteOperator does not correctly report its inversion options!");
      if (opts.template get<std::string>("type") == type)
        return opts;
    }
    DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong,
               "Given type is not one of invert_options()!"
                   << "\n"
                   << "   type = " << type << "\n"
                   << "   invert_options() = " << print(this->invert_options()));
    return XT::Common::Configuration();
  } // ... invert_options(...)

protected:
  /**
   * Used by derived DiscreteOperators to feed the implementations of invert_options().
   * Simply return a vector of configurations with each at least a "type" key.
   */
  virtual std::vector<XT::Common::Configuration> all_invert_options() const
  {
    auto newton_opts = default_newton_solve_options();
    newton_opts["type"] = "newton";
    return {newton_opts};
  }

public:
  /**
   * \brief Computes the inverse action of the operator.
   *
   * Currently implemented is the dampened Newton from [DF2015, Sec. 8.4.4.1].
   *
   * \note Given source is used as initial guess.
   *
   * \todo Allow to pass jacobian options as subcfg in newton.
   **/
  virtual void apply_inverse(const VectorType& range_vector,
                             VectorType& source_vector,
                             const XT::Common::Configuration& opts,
                             const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply_inverse(range_vector.sup_norm()=" << range_vector.sup_norm()
                << ", source_vector.sup_norm()=" << source_vector.sup_norm()
                << ",\n   opts=" << print(opts, {{"oneline", "true"}}) << ",\n   param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    this->assert_matching_source(source_vector);
    this->assert_apply_inverse_opts(opts);
    const auto type = opts.get<std::string>("type");
    if (type == "newton") {
      LOG_(info) << "computing inverse by dampened newton ..." << std::endl;
      newton_solve(*this, range_vector, source_vector, param, opts);
    } else
      DUNE_THROW(Exceptions::operator_error,
                 "unknown apply_inverse type requested (" << type << "),"
                                                          << "\n"
                                                          << "   available types are " << print(this->invert_options())
                                                          << std::endl);
  } // ... apply_inverse(...)

  /// \}
  /// \name These apply methods are provided for convenience (discrete source variants).
  /// \{

  virtual VectorType apply(const VectorType& source_vector, const XT::Common::Parameter& param = {}) const
  {
    VectorType range(this->range_space().mapper().size(), 0.);
    this->apply(source_vector, range, param);
    return range;
  }

  virtual void apply(const ConstDiscreteSourceFunctionType& discrete_source_function,
                     DiscreteRangeFunctionType& discrete_range_function,
                     const XT::Common::Parameter& param = {}) const
  {
    this->apply(discrete_source_function.dofs().vector(), discrete_range_function.dofs().vector(), param);
  }

  virtual DiscreteRangeFunctionType apply(const ConstDiscreteSourceFunctionType& discrete_source_function,
                                          const XT::Common::Parameter& param = {}) const
  {
    DiscreteRangeFunctionType discrete_range_function(this->range_space());
    this->apply(discrete_source_function, discrete_range_function, param);
    return discrete_range_function;
  }

  /// \}
  /// \name These jacobian variants are provided for convenience (vector variants, simplified type/options).
  /// \{

  virtual void jacobian(const VectorType& source_vector,
                        MatrixOperatorType& jacobian_op,
                        const std::string& type,
                        const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm() << ",\n   type=" << type
                << ", param=" << print(param) << std::endl;
    return this->jacobian(source_vector, jacobian_op, this->jacobian_options(type), param);
  }

  virtual void jacobian(const VectorType& source_vector,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm() << ", param=" << print(param)
                << std::endl;
    return this->jacobian(source_vector, jacobian_op, this->jacobian_options().at(0), param);
  }

  /// \}
  /// \name These jacobian variants are provided for convenience (vector variants, creates matching return op).
  /// \{

  virtual MatrixOperatorType jacobian(const VectorType& source_vector,
                                      const XT::Common::Configuration& opts,
                                      const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << print(param) << std::endl;
    auto jacobian_op = this->empty_jacobian_op();
    this->jacobian(source_vector, jacobian_op, opts, param);
    return jacobian_op;
  }

  virtual MatrixOperatorType
  jacobian(const VectorType& source_vector, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm() << ", type=" << type
                << ", param=" << print(param) << std::endl;
    auto jacobian_op = this->empty_jacobian_op();
    this->jacobian(source_vector, jacobian_op, type, param);
    return jacobian_op;
  }

  virtual MatrixOperatorType jacobian(const VectorType& source_vector, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm() << ", param=" << print(param)
                << std::endl;
    auto jacobian_op = this->empty_jacobian_op();
    this->jacobian(source_vector, jacobian_op, param);
    return jacobian_op;
  }

  /// \}
  /// \name These jacobian variants are provided for convenience (discrete function variants).
  /// \{

  virtual void jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Configuration& opts,
                        const XT::Common::Parameter& param = {}) const
  {
    this->jacobian(discrete_source_function.dofs().vector(), jacobian_op, opts, param);
  }

  virtual void jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                        MatrixOperatorType& jacobian_op,
                        const std::string& type,
                        const XT::Common::Parameter& param = {}) const
  {
    this->jacobian(discrete_source_function.dofs().vector(), jacobian_op, type, param);
  }

  virtual void jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Parameter& param = {}) const
  {
    this->jacobian(discrete_source_function.dofs().vector(), jacobian_op, param);
  }

  virtual MatrixOperatorType jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                                      const XT::Common::Configuration& opts,
                                      const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(discrete_source_function.dofs().vector(), opts, param);
  }

  virtual MatrixOperatorType jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                                      const std::string& type,
                                      const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(discrete_source_function.dofs().vector(), type, param);
  }

  virtual MatrixOperatorType jacobian(const ConstDiscreteSourceFunctionType& discrete_source_function,
                                      const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(discrete_source_function.dofs().vector(), param);
  }

  /// \}
  /// \name These apply_inverse variants are provided for convenience (vector variants, simplified type/options).
  /// \{

  virtual void apply_inverse(const VectorType& range_vector,
                             VectorType& source_vector,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range_vector, source_vector, this->invert_options(type), param);
  }

  virtual void apply_inverse(const VectorType& range_vector,
                             VectorType& source_vector,
                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(this->invert_options().empty(),
                  Exceptions::operator_error,
                  "invert_options must not be empty to support apply_inverse!");
    this->apply_inverse(range_vector, source_vector, this->invert_options().at(0), param);
  }

  /// \}
  /// \name These apply_inverse variants are provided for convenience (vector variants, creates matching return vector).
  /// \{

  virtual VectorType apply_inverse(const VectorType& range_vector,
                                   const XT::Common::Configuration& opts,
                                   const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range_vector, source, opts, param);
    return source;
  }

  virtual VectorType
  apply_inverse(const VectorType& range_vector, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range_vector, source, type, param);
    return source;
  }

  virtual VectorType apply_inverse(const VectorType& range_vector, const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range_vector, source, param);
    return source;
  }

  /// \}
  /// \name These apply_inverse variants are provided for convenience (discrete function variants).
  /// \{

  virtual void apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                             DiscreteSourceFunctionType& discrete_source_function,
                             const XT::Common::Configuration& opts,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(discrete_range_function.dofs().vector(), discrete_source_function.dofs().vector(), opts, param);
  }

  virtual void apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                             DiscreteSourceFunctionType& discrete_source_function,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(discrete_range_function.dofs().vector(), discrete_source_function.dofs().vector(), type, param);
  }

  virtual void apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                             DiscreteSourceFunctionType& discrete_source_function,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(discrete_range_function.dofs().vector(), discrete_source_function.dofs().vector(), param);
  }

  virtual DiscreteSourceFunctionType apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                                                   const XT::Common::Configuration& opts,
                                                   const XT::Common::Parameter& param = {}) const
  {
    return DiscreteSourceFunctionType(this->source_space(),
                                      this->apply_inverse(discrete_range_function.dofs().vector(), opts, param));
  }

  virtual DiscreteSourceFunctionType apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                                                   const std::string& type,
                                                   const XT::Common::Parameter& param = {}) const
  {
    return DiscreteSourceFunctionType(this->source_space(),
                                      this->apply_inverse(discrete_range_function.dofs().vector(), type, param));
  }

  virtual DiscreteSourceFunctionType apply_inverse(const ConstDiscreteRangeFunctionType& discrete_range_function,
                                                   const XT::Common::Parameter& param = {}) const
  {
    return DiscreteSourceFunctionType(this->source_space(),
                                      this->apply_inverse(discrete_range_function.dofs().vector(), param));
  }

  /// \}
  /// \name const operator* and operator/ variants
  /// \{

  virtual ConstLincombOperatorType operator*(const FieldType& alpha) const
  {
    return make_operator_mul(*this, alpha);
  }

  virtual ConstLincombOperatorType operator/(const FieldType& alpha) const
  {
    return make_operator_div(*this, alpha);
  }

  /// \}
  /// \name mutable operator* and operator/ variants
  /// \{

  virtual LincombOperatorType operator*(const FieldType& alpha)
  {
    return make_operator_mul(*this, alpha);
  }

  virtual LincombOperatorType operator/(const FieldType& alpha)
  {
    return make_operator_div(*this, alpha);
  }

  /// \{
  /// \name const operator+ variants
  /// \{

  virtual ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const
  {
    return make_operator_addsub(*this, other, /*add=*/true);
  }

  virtual ConstLincombOperatorType operator+(const ThisType& other) const
  {
    return make_operator_addsub(*this, other, /*add=*/true);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  virtual ConstLincombOperatorType operator+(const VectorType& vector) const
  {
    return make_operator_addsub(*this, vector, /*add=*/true);
  }

  /// \}
  /// \name mutable operator+ variants
  /// \{

  virtual LincombOperatorType operator+(LincombOperatorType& other)
  {
    return make_operator_addsub(*this, other, /*add=*/true);
  }

  virtual LincombOperatorType operator+(ThisType& other)
  {
    return make_operator_addsub(*this, other, /*add=*/true);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  virtual LincombOperatorType operator+(const VectorType& vector)
  {
    return make_operator_addsub(*this, vector, /*add=*/true);
  }

  /// \}
  /// \name const operator- variants
  /// \{

  virtual ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const
  {
    return make_operator_addsub(*this, other, /*add=*/false);
  }

  virtual ConstLincombOperatorType operator-(const ThisType& other) const
  {
    return make_operator_addsub(*this, other, /*add=*/false);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  virtual ConstLincombOperatorType operator-(const VectorType& vector) const
  {
    return make_operator_addsub(*this, vector, /*add=*/false);
  }

  /// \}
  /// \name mutable operator- variants
  /// \{

  virtual LincombOperatorType operator-(LincombOperatorType& other)
  {
    return make_operator_addsub(*this, other, /*add=*/false);
  }

  virtual LincombOperatorType operator-(ThisType& other)
  {
    return make_operator_addsub(*this, other, /*add=*/false);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  virtual LincombOperatorType operator-(const VectorType& vector)
  {
    return make_operator_addsub(*this, vector, /*add=*/false);
  }

  /// \}

  MatrixOperatorType empty_jacobian_op() const
  {
    std::string prefix = this->logger.prefix + "_jac";
    for (auto&& character : {"(", ")", "+", "-", "*", "/"})
      if (this->logger.prefix.find(character) != std::string::npos)
        prefix = "(" + this->logger.prefix + ")_jac";
    return MatrixOperatorType(
        this->assembly_grid_view(),
        this->source_space(),
        this->range_space(),
        new MatrixType(this->range_space().mapper().size(),
                       this->source_space().mapper().size(),
                       make_sparsity_pattern(this->range_space(), this->source_space(), this->assembly_grid_view())),
        prefix,
        this->logger.state);
  } // ... empty_jacobian_op(...)

protected:
  void assert_matching_source(const VectorType& source_vector) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source_vector),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "source_space().mapper().size() = " << this->source_space().mapper().size() << "\n"
                                                      << "   source_vector.size() = " << source_vector.size());
  } // ... assert_matching_source(...)

  void assert_matching_source(const DiscreteSourceFunctionType& discrete_source_function) const
  {
    DUNE_THROW_IF(!this->source_space().contains(discrete_source_function),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "source_space().mapper().size() = "
                      << this->source_space().mapper().size() << "\n"
                      << "   discrete_source_function.dofs().vector().size() = "
                      << discrete_source_function.dofs().vector().size() << "\n"
                      << "   source_space.type() = " << this->source_space().type() << "\n"
                      << "   discrete_source_function.space().type() = " << discrete_source_function.space().type());
  } // ... assert_matching_source(...)

  void assert_jacobian_opts(const XT::Common::Configuration& opts) const
  {
    const auto available_types = this->jacobian_options();
    DUNE_THROW_IF(!opts.has_key("type"),
                  Exceptions::operator_error,
                  "missing key 'type' in given opts!"
                      << "\n\n   available type: " << print(available_types) << "\n   opts = " << print(opts));
    const auto type = opts.get<std::string>("type");
    DUNE_THROW_IF(std::find(available_types.begin(), available_types.end(), type) == available_types.end(),
                  Exceptions::operator_error,
                  "requested jacobian type is not one of the available ones!"
                      << "\n\n   type = " << type << "\n   jacobian_options() = " << print(available_types));
  } // ... assert_jacobian_opts(...)

  void assert_apply_inverse_opts(const XT::Common::Configuration& opts) const
  {
    const auto available_types = this->invert_options();
    DUNE_THROW_IF(!opts.has_key("type"),
                  Exceptions::operator_error,
                  "missing key 'type' in given opts!"
                      << "\n\n   available type: " << print(available_types) << "\n   opts = " << print(opts));
    const auto type = opts.get<std::string>("type");
    DUNE_THROW_IF(std::find(available_types.begin(), available_types.end(), type) == available_types.end(),
                  Exceptions::operator_error,
                  "requested inversion type is not one of the available ones!"
                      << "\n\n   type = " << type << "\n   invert_options() = " << print(available_types));
  } // ... assert_jacobian_opts(...)

  // used in derived classes, template required for add(self, ...) to resolve correctly
  template <typename SelfType>
  static auto make_operator_mul(SelfType& self, const FieldType& alpha)
  {
    constexpr bool is_mutable = std::is_same_v<SelfType, std::remove_const_t<SelfType>>;
    using ReturnType = std::conditional_t<is_mutable, LincombOperatorType, ConstLincombOperatorType>;
    if (self.logger.debug_enabled()) // cannot use LOG_(debug) here
      self.logger.debug() << "operator*(this=" << (is_mutable ? "mutable" : "const") << ", alpha=" << alpha << ")"
                          << std::endl;
    ReturnType ret(self.assembly_grid_view(),
                   self.source_space(),
                   self.range_space(),
                   "(" + self.logger.prefix + ")*" + XT::Common::to_string(alpha),
                   self.logger.state);
    ret.add(self, alpha);
    return ret;
  } // ... make_operator_mul(...)

  // used in derived classes, template required for add(self, ...) to resolve correctly
  template <typename SelfType>
  static auto make_operator_div(SelfType& self, const FieldType& alpha)
  {
    constexpr bool is_mutable = std::is_same_v<SelfType, std::remove_const_t<SelfType>>;
    using ReturnType = std::conditional_t<is_mutable, LincombOperatorType, ConstLincombOperatorType>;
    if (self.logger.debug_enabled()) // cannot use LOG_(debug) here
      self.logger.debug() << "operator/(this=" << (is_mutable ? "mutable" : "const") << ", alpha=" << alpha << ")"
                          << std::endl;
    ReturnType ret(self.assembly_grid_view(),
                   self.source_space(),
                   self.range_space(),
                   "(" + self.logger.prefix + ")/" + XT::Common::to_string(alpha),
                   self.logger.state);
    ret.add(self, 1. / alpha);
    return ret;
  } // ... make_operator_div(...)

  // used in derived classes, template required for add(self, ...) to resolve correctly
  template <typename SelfType, typename OtherType>
  static auto make_operator_addsub(SelfType& self, OtherType& other, const bool add)
  {
    constexpr bool is_mutable = std::is_same_v<SelfType, std::remove_const_t<SelfType>>;
    using ReturnType = std::conditional_t<is_mutable, LincombOperatorType, ConstLincombOperatorType>;
    constexpr bool other_is_lincomb = std::is_same_v<std::remove_const_t<OtherType>, ReturnType>;
    if (self.logger.debug_enabled()) // cannot use LOG_(debug) here
      self.logger.debug() << "operator" << (add ? "+" : "-") << "(this=" << (is_mutable ? "mutable" : "const")
                          << ", other_" << (is_mutable ? "" : "const_") << (other_is_lincomb ? "lincomb_" : "")
                          << "op=" << &other << ")" << std::endl;
    ReturnType ret(self.assembly_grid_view(),
                   self.source_space(),
                   self.range_space(),
                   self.logger.prefix + " " + (add ? "+" : "-") + " " + other.logger.prefix,
                   self.logger.get_state_or(other.logger.state));
    ret.add(self, 1.);
    ret.add(other, add ? 1. : -1.);
    return ret;
  } // ... make_operator_addsub(...)

  // used in derived classes, template required for add(self, ...) to resolve correctly
  template <typename SelfType>
  static auto make_operator_addsub(SelfType& self, const VectorType& vector, const bool add)
  {
    constexpr bool is_mutable = std::is_same_v<SelfType, std::remove_const_t<SelfType>>;
    using ReturnType = std::conditional_t<is_mutable, LincombOperatorType, ConstLincombOperatorType>;
    if (self.logger.debug_enabled()) // cannot use LOG_(debug) here
      self.logger.debug() << "operator" << (add ? "+" : "-") << "(this=" << (is_mutable ? "mutable" : "const")
                          << ", vector.sup_norm()=" << vector.sup_norm() << ")" << std::endl;
    const std::string id = XT::Common::is_zero(vector.sup_norm()) ? "0" : "vec";
    ReturnType ret(self.assembly_grid_view(),
                   self.source_space(),
                   self.range_space(),
                   self.logger.prefix + " " + (add ? "+" : "-") + " " + id,
                   self.logger.state);
    ret.add(self, 1.);
    ret.add(new ConstantOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>(
                self.assembly_grid_view(), self.source_space(), self.range_space(), vector, id, self.logger.state),
            add ? 1. : -1.);
    return ret;
  } // ... make_operator_addsub(...)
}; // class OperatorInterface


} // namespace GDT
} // namespace Dune

#include "constant.hh"
#include "lincomb.hh"
#include "matrix.hh"

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
