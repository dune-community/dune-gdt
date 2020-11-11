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

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/interpolations.hh>
#include <dune/gdt/print.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

namespace Dune {
namespace GDT {


//// forwards, required for operator +-*/
// template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
// class ConstantOperator;

// template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
// class ConstLincombOperator;

// template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
// class LincombOperator;


// forward, required for the jacobian
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
  static_assert(XT::Grid::is_view<SourceGridView>::value, "");
  static_assert(XT::Grid::is_view<RangeGridView>::value, "");
  static_assert(SourceGridView::dimension == RangeGridView::dimension, "");

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
                                 const std::array<bool, 3>& logging_enabled = {false, false, true})
    : XT::Common::ParametricInterface(param_type)
    , Logger(logging_prefix.empty() ? "BilinearFormInterface" : logging_prefix, logging_enabled)
  {
    LOG_(info) << "BilinearFormInterface(param_type=" << param_type << ")" << std::endl;
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
  virtual ThisType& assemble(const bool use_tbb = false)
  {
    LOG_(debug) << "assemble(use_tbb=" << use_tbb << ")" << std::endl;
    LOG_(info) << "not assembling bilinear form" << std::endl;
    return *this;
  }

  FieldType norm(RangeFunctionType range_function, const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "norm(range_function=" << &range_function << ", param=" << param << ")" << std::endl;
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
class OperatorInterface : public BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>
{
  static_assert(XT::LA::is_vector<Vector>::value, "");

  using ThisType = OperatorInterface;
  using BaseType = BilinearFormInterface<SGV, s_r, s_rC, r_r, r_rC, F, RGV>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::RangeFunctionType;
  using typename BaseType::SourceFunctionType;

  using VectorType = Vector;
  using V = VectorType;

  using RangeSpaceType = SpaceInterface<RGV, r_r, r_rC, F>;
  using DiscreteRangeFunctionType = DiscreteFunction<V, RGV, r_r, r_rC, F>;

  explicit OperatorInterface(const XT::Common::ParameterType& param_type = {},
                             const std::string& logging_prefix = "",
                             const bool logging_disabled = true)
    : BaseType(param_type, logging_prefix.empty() ? "Operator" : logging_prefix, logging_disabled)
  {
    LOG_(info) << "Operator(param_type=" << param_type << ")" << std::endl;
  }

  OperatorInterface(const ThisType& other) = default;

  OperatorInterface(ThisType&& source) = default;

  virtual ~OperatorInterface() = default;

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
                << ", param=" << param << ")" << std::endl;
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
                << ", source_function=" << &source_function << ", param=" << param << ")" << std::endl;
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
                << ", discrete_range_function=" << &discrete_range_function << ", param=" << param << ")"
                << "\n"
                << "redicrecting to discrete apply() variant ..." << std::endl;
    this->apply(source_function, discrete_range_function.dofs().vector(), param);
  }

  virtual DiscreteRangeFunctionType apply(SourceFunctionType source_function,
                                          const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "apply(source_function=" << &source_function << ", param=" << param << ")"
                << "\n"
                << "creating discrete_range_function and redicrecting to discrete apply() variant ..." << std::endl;
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
}; // class OperatorInterface


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
class DiscreteOperatorInterface : public OperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, XT::LA::vector_t<Matrix>, RGV>
{
  static_assert(XT::Grid::is_view<AssemblyGridView>::value, "");
  static_assert(XT::LA::is_matrix<Matrix>::value, "");

  using ThisType = DiscreteOperatorInterface;
  using BaseType = OperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, XT::LA::vector_t<Matrix>, RGV>;

public:
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

  explicit DiscreteOperatorInterface(const XT::Common::ParameterType& param_type = {},
                                     const std::string& logging_prefix = "",
                                     const bool logging_disabled = true)
    : BaseType(param_type, logging_prefix.empty() ? "DiscreteOperator" : logging_prefix, logging_disabled)
  {
    LOG_(info) << "DiscreteOperator(param_type=" << param_type << ")" << std::endl;
  }

  DiscreteOperatorInterface(const ThisType& other) = default;

  DiscreteOperatorInterface(ThisType&& source) = default;

  virtual ~DiscreteOperatorInterface() = default;

  // pull in methods from BilinearFormInterface and OperatorInterface
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
  /// \name These methods are required by OperatorInterface and default-implemented here.
  /// \{

  /// \brief Carries out an interpolation of source_function and calls the other apply() afterwards (not optimal!).
  virtual void apply(SourceFunctionType source_function,
                     VectorType& range_vector,
                     const XT::Common::Parameter& param = {}) const override
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << param << ")" << std::endl;
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
                << ", source_vector.sup_norm()=" << source_vector.sup_norm() << ", param=" << param << ")" << std::endl;
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
    LOG_(debug) << "norm(range_vector.sup_norm=" << range_vector.sup_norm() << ", param=" << param << ")" << std::endl;
    if constexpr ((r_r == s_r) && (r_rC == s_rC)) {
      this->assert_matching_range(range_vector);
      this->assert_matching_source(range_vector);
      LOG_(debug) << "calling norm(discrete_range_function, param)... " << std::endl;
      // sanity checks for zero are carried out in the BilinearForm
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
                << ",\n   opts=" << print(opts, {{"oneline", "true"}}) << ",\n   param=" << param << ")" << std::endl;
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
    return {{{"type", "newton"}, {"precision", "1e-7"}, {"max_iter", "100"}, {"max_dampening_iter", "1000"}}};
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
#if 0
    LOG_(debug) << "apply_inverse(range_vector.sup_norm()=" << range_vector.sup_norm()
                << ", source_vector.sup_norm()=" << source_vector.sup_norm()
                << ",\n   opts=" << print(opts, {{"oneline", "true"}}) << ",\n   param=" << param << ")" << std::endl;
    this->assert_matching_range(range_vector);
    this->assert_matching_source(source_vector);
    this->assert_apply_inverse_opts(opts);
    const auto type = opts.get<std::string>("type");
    const XT::Common::Configuration default_opts = this->invert_options(type);
    if (type == "newton") {
      // some preparations
      auto residual_op = *this - range_vector;
      auto residual = range_vector.copy();
      auto update = source_vector.copy();
      auto candidate = source_vector.copy();
      // one matrix for all jacobians
      MatrixOperatorType jacobian_op(
            this->assembly_grid_view(),
            this->source_space(),
            this->range_space(),
            new MatrixType(this->range_space().mapper().size(),
                           this->source_space().mapper().size(),
                           make_element_and_intersection_sparsity_pattern(
                               this->range_space(), this->source_space(), this->source_space().grid_view())));
      XT::LA::Solver<M> jacobian_solver(jacobian_op.matrix());
      const auto precision = opts.get("precision", default_opts.get<double>("precision"));
      const auto max_iter = opts.get("max_iter", default_opts.get<size_t>("max_iter"));
      const auto max_dampening_iter = opts.get("max_dampening_iter", default_opts.get<size_t>("max_iter"));
      LOG_(info) << "computing inverse by dampened newton ..." << std::endl;
      size_t l = 0;
      Timer timer;
      while (true) {
        timer.reset();
        LOG_(debug) << "l = " << l << ": computing residual ... " << std::flush;
        residual_op.apply(source_vector, residual, param);
        auto res = residual.l2_norm();
        LOG_(debug) << "took " << timer.elapsed() << "s, |residual|_l2 = " << res << std::endl;
        if (res < precision) {
          LOG_(debug) << "       residual below tolerance, succeeded!" << std::endl;
          break;
        }
        DUNE_THROW_IF(l >= max_iter,
                      Exceptions::operator_error,
                      "max iterations reached!\n|residual|_l2 = " << res << "\nopts:\n"
                                                                  << opts);
        LOG_(debug) << "       computing jacobi matrix ... " << std::flush;
        timer.reset();
        jacobian_op.matrix() *= 0.;
        residual_op.jacobian(source_vector, jacobian_op, {{"type", residual_op.jacobian_options().at(0)}}, param);
        jacobian_op.assemble(/*use_tbb=*/true);
        LOG_(debug) << "took " << timer.elapsed() << "s\n       solving for defect ... " << std::flush;
        timer.reset();
        residual *= -1.;
        update = source_vector; // <- initial guess for the linear solver
        bool linear_solve_succeeded = false;
        std::vector<std::string> tried_linear_solvers;
        for (const auto& linear_solver_type : jacobian_solver.types()) {
          try {
            tried_linear_solvers.push_back(linear_solver_type);
            jacobian_solver.apply(
                residual,
                update,
                {{"type", linear_solver_type}, {"precision", XT::Common::to_string(0.1 * precision)}});
            linear_solve_succeeded = true;
          } catch (const XT::LA::Exceptions::linear_solver_failed&) {
          }
        }
        DUNE_THROW_IF(!linear_solve_succeeded,
                      Exceptions::operator_error,
                      "could not solve linear system for defect!\nTried the following linear solvers: "
                          << tried_linear_solvers << "\nopts:\n"
                          << opts);
        LOG_(debug) << "took " << timer.elapsed() << "s";
        if (tried_linear_solvers.size() > 1) {
          LOG_(debug) << ", took " << tried_linear_solvers.size() << " attempts with different linear solvers";
        }
        LOG_(debug) << "\n       computing update ... " << std::flush;
        timer.reset();
        // try the automatic dampening strategy proposed in [DF2015, Sec. 8.4.4.1, p. 432]
        size_t k = 0;
        auto candidate_res = 2 * res; // any number such that we enter the while loop at least once
        double lambda = 1;
        while (!(candidate_res / res < 1)) {
          DUNE_THROW_IF(k >= max_dampening_iter,
                        Exceptions::operator_error,
                        "max iterations reached when trying to compute automatic dampening!\n|residual|_l2 = "
                            << res << "\nl = " << l << "\nopts:\n"
                            << opts);
          candidate = source_vector + update * lambda;
          residual_op.apply(candidate, residual, param);
          candidate_res = residual.l2_norm();
          lambda /= 2;
          k += 1;
        }
        source_vector = candidate;
        LOG_(debug) << "took " << timer.elapsed() << "s and a dampening of " << 2 * lambda << std::endl;
        l += 1;
      }
    } else
      DUNE_THROW(Exceptions::operator_error,
                 "unknown apply_inverse type requested (" << type << "),"
                                                          << "\n"
                                                          << "   available types are " << print(this->invert_options())
                                                          << std::endl);
#endif // 0
  } // ... apply_inverse(...)

  /// \}
  /// \name These apply methods are provided for convenience (non-discrete source function variants).
  /// \{

  virtual void apply(SourceFunctionType source_function,
                     DiscreteRangeFunctionType& discrete_range_function,
                     const XT::Common::Parameter& param = {}) const
  {
    this->apply(source_function, discrete_range_function.dofs().vector(), param);
  }

  virtual DiscreteRangeFunctionType apply(SourceFunctionType source_function,
                                          const XT::Common::Parameter& param = {}) const
  {
    DiscreteRangeFunctionType discrete_range_function(this->range_space());
    this->apply(source_function, discrete_range_function.dofs().vector(), param);
    return discrete_range_function;
  }

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
                << ", param=" << param << std::endl;
    return this->jacobian(source_vector, jacobian_op, this->jacobian_options(type), param);
  }

  virtual void jacobian(const VectorType& source_vector,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Parameter& param = {}) const
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm() << ", param=" << param
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
    std::string derived_logging_prefix = "";
    if (this->logger.info_enabled) {
      derived_logging_prefix = this->logger.prefix + "_jac";
      this->logger.debug() << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                           << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << param << std::endl;
    }
    auto jacobian_op = this->empty_jacobian_op(derived_logging_prefix);
    this->jacobian(source_vector, jacobian_op, opts, param);
    return jacobian_op;
  } // ... jacobian(...)

  virtual MatrixOperatorType
  jacobian(const VectorType& source_vector, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.info_enabled) {
      derived_logging_prefix = this->logger.prefix + "_jac";
      this->logger.info() << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm() << ", type=" << type
                          << ", param=" << param << std::endl;
    }
    auto jacobian_op = this->empty_jacobian_op(derived_logging_prefix);
    this->jacobian(source_vector, jacobian_op, type, param);
    return jacobian_op;
  } // ... jacobian(...)

  virtual MatrixOperatorType jacobian(const VectorType& source_vector, const XT::Common::Parameter& param = {}) const
  {
    std::string derived_logging_prefix = "";
    if (this->logger.info_enabled) {
      derived_logging_prefix = this->logger.prefix + "_jac";
      this->logger.info() << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm() << ", param=" << param
                          << std::endl;
    }
    auto jacobian_op = this->empty_jacobian_op(derived_logging_prefix);
    this->jacobian(source_vector, jacobian_op, param);
    return jacobian_op;
  } // ... jacobian(...)

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

protected:
  MatrixOperatorType empty_jacobian_op(const std::string& logging_prefix = "") const
  {
    return MatrixOperatorType(
        this->assembly_grid_view(),
        this->source_space(),
        this->range_space(),
        new MatrixType(this->range_space().mapper().size(),
                       this->source_space().mapper().size(),
                       make_element_and_intersection_sparsity_pattern(
                           this->range_space(), this->source_space(), this->assembly_grid_view())),
        logging_prefix);
  } // ... empty_jacobian_op(...)

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
}; // class DiscreteOperatorInterface


//  virtual ConstLincombOperatorType operator*(const FieldType& alpha) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = this->logger.prefix + "*" + XT::Common::to_string(alpha);
//      this->logger.debug() << "operator*(alpha=" << alpha << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, alpha);
//    return ret;
//  }

//  virtual LincombOperatorType operator*(const FieldType& alpha)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = this->logger.prefix + "*" + XT::Common::to_string(alpha);
//      this->logger.debug() << "operator*(alpha=" << alpha << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, alpha);
//    return ret;
//  }

//  virtual ConstLincombOperatorType operator/(const FieldType& alpha) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = this->logger.prefix + "/" + XT::Common::to_string(alpha);
//      this->logger.debug() << "operator/(alpha=" << alpha << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1. / alpha);
//    return ret;
//  }

//  virtual LincombOperatorType operator/(const FieldType& alpha)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = this->logger.prefix + "/" + XT::Common::to_string(alpha);
//      this->logger.debug() << "operator/(alpha=" << alpha << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1. / alpha);
//    return ret;
//  }

//  /// \name const operator+ variants
//  /// \{

//  virtual ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "ConstLincombOperator";
//      this->logger.debug() << "operator+(other_const_lincomb_op=" << &other << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, 1.);
//    return ret;
//  }

//  virtual ConstLincombOperatorType operator+(const ThisType& other) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "ConstLincombOperator";
//      this->logger.debug() << "operator+(other_const_op=" << &other << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, 1.);
//    return ret;
//  }

//  /**
//   * \note vector is interpreted as a ConstantOperator
//   * \sa ConstantOperator
//   */
//  virtual ConstLincombOperatorType operator+(const VectorType& vector) const
//  {
//    std::string derived_logging_prefix_clop = "";
//    std::string derived_logging_prefix_cop = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix_clop = "ConstLincombOperator";
//      derived_logging_prefix_cop = "ConstantOperator";
//      this->logger.debug() << "operator+(vector.sup_norm()=" << vector.sup_norm() << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix_clop);
//    ret.add(*this, 1.);
//    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(
//                this->source_space(), this->range_space(), vector, derived_logging_prefix_cop),
//            1.);
//    return ret;
//  }

//  /// \}
//  /// \name mutable operator+ variants
//  /// \{

//  virtual LincombOperatorType operator+(LincombOperatorType& other)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "LincombOperator";
//      this->logger.debug() << "operator+(other_lincomb_op=" << &other << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, 1.);
//    return ret;
//  }

//  virtual LincombOperatorType operator+(ThisType& other)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "LincombOperator";
//      this->logger.debug() << "operator+(other_op=" << &other << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, 1.);
//    return ret;
//  }

//  /**
//   * \note vector is interpreted as a ConstantOperator
//   * \sa ConstantOperator
//   */
//  virtual LincombOperatorType operator+(const VectorType& vector)
//  {
//    std::string derived_logging_prefix_lop = "";
//    std::string derived_logging_prefix_cop = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix_lop = "LincombOperator";
//      derived_logging_prefix_cop = "ConstantOperator";
//      this->logger.debug() << "operator+(vector.sup_norm()=" << vector.sup_norm() << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix_lop);
//    ret.add(*this, 1.);
//    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(
//                this->source_space(), this->range_space(), vector, derived_logging_prefix_cop),
//            1.);
//    return ret;
//  }

//  /// \}
//  /// \name const operator- variants
//  /// \{

//  virtual ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "ConstLincombOperator";
//      this->logger.debug() << "operator-(other_const_lincomb_op=" << &other << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, -1.);
//    return ret;
//  }

//  virtual ConstLincombOperatorType operator-(const ThisType& other) const
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "ConstLincombOperator";
//      this->logger.debug() << "operator-(other_op=" << &other << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, -1.);
//    return ret;
//  }

//  /**
//   * \note vector is interpreted as a ConstantOperator
//   * \sa ConstantOperator
//   */
//  virtual ConstLincombOperatorType operator-(const VectorType& vector) const
//  {
//    std::string derived_logging_prefix_clop = "";
//    std::string derived_logging_prefix_cop = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix_clop = "ConstLincombOperator";
//      derived_logging_prefix_cop = "ConstantOperator";
//      this->logger.debug() << "operator-(vector.sup_norm()=" << vector.sup_norm() << ")" << std::endl;
//    }
//    ConstLincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix_clop);
//    ret.add(*this, 1.);
//    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(
//                this->source_space(), this->range_space(), vector, derived_logging_prefix_cop),
//            -1.);
//    return ret;
//  }

//  /// \}
//  /// \name mutable operator- variants between arbitrary operators
//  /// \{

//  virtual LincombOperatorType operator-(LincombOperatorType& other)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "LincombOperator";
//      this->logger.debug() << "operator-(other_lincomb_op=" << &other << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, -1.);
//    return ret;
//  }

//  virtual LincombOperatorType operator-(ThisType& other)
//  {
//    std::string derived_logging_prefix = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix = "LincombOperator";
//      this->logger.debug() << "operator-(other_op=" << &other << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix);
//    ret.add(*this, 1.);
//    ret.add(other, -1.);
//    return ret;
//  }

//  /**
//   * \note vector is interpreted as a ConstantOperator
//   * \sa ConstantOperator
//   */
//  virtual LincombOperatorType operator-(const VectorType& vector)
//  {
//    std::string derived_logging_prefix_lop = "";
//    std::string derived_logging_prefix_cop = "";
//    if (this->logger.debug_enabled) {
//      derived_logging_prefix_lop = "LincombOperator";
//      derived_logging_prefix_cop = "ConstantOperator";
//      this->logger.debug() << "operator-(vector.sup_norm()=" << vector.sup_norm() << ")" << std::endl;
//    }
//    LincombOperatorType ret(this->source_space(), this->range_space(), derived_logging_prefix_lop);
//    ret.add(*this, 1.);
//    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(
//                this->source_space(), this->range_space(), vector, derived_logging_prefix_cop),
//            -1.);
//    return ret;
//  }

/// \}


} // namespace GDT
} // namespace Dune

#include "constant.hh"
//#include "lincomb.hh"
#include "matrix-based.hh"

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
