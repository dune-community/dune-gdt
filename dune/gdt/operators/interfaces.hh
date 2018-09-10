// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reseVed.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016, 2018)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_OPERATORS_INTERFACES_HH
#define DUNE_GDT_OPERATORS_INTERFACES_HH

#include <cmath>
#include <type_traits>

#include <dune/xt/common/parameter.hh>
#include <dune/xt/common/configuration.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/tools/sparsity-pattern.hh>

namespace Dune {
namespace GDT {
namespace internal {


template <class Matrix,
          class SourceGridView,
          size_t source_dim,
          size_t source_dim_cols,
          size_t range_dim,
          size_t range_dim_cols,
          class RangeGridView>
class AssertArgumentsOfOperatorInterface
{
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_view<SourceGridView>::value, "");
  static_assert(source_dim > 0, "");
  static_assert(source_dim_cols > 0, "");
  static_assert(range_dim > 0, "");
  static_assert(range_dim_cols > 0, "");
  static_assert(XT::Grid::is_view<RangeGridView>::value, "");
  static_assert(SourceGridView::dimension == RangeGridView::dimension, "");
}; // class AssertArgumentsOfOperatorInterface


} // namespace internal


// forwards, required for operator +-*/
template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
class ConstLincombOperator;

template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
class LincombOperator;


// forward, required for the jacobian
template <class M, class SGV, size_t s_r, size_t s_rC, size_t r_r, size_t r_rC, class RGV>
class MatrixOperator;


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
template <class Matrix,
          class SourceGridView,
          size_t source_dim = 1,
          size_t source_dim_cols = 1,
          size_t range_dim = source_dim,
          size_t range_dim_cols = source_dim_cols,
          class RangeGridView = SourceGridView>
class OperatorInterface : internal::AssertArgumentsOfOperatorInterface<Matrix,
                                                                       SourceGridView,
                                                                       source_dim,
                                                                       source_dim_cols,
                                                                       range_dim,
                                                                       range_dim_cols,
                                                                       RangeGridView>,
                          public XT::Common::ParametricInterface
{
public:
  using MatrixType = Matrix;
  using VectorType = XT::LA::vector_t<MatrixType>;
  using FieldType = typename MatrixType::ScalarType;

  using M = Matrix;
  using V = VectorType;
  using F = FieldType;

  using SGV = SourceGridView;
  static const constexpr size_t s_r = source_dim;
  static const constexpr size_t s_rC = source_dim_cols;

  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_dim;
  static const constexpr size_t r_rC = range_dim_cols;

  using SourceSpaceType = SpaceInterface<SGV, s_r, s_rC, F>;
  using SourceFunctionType = DiscreteFunction<V, SGV, s_r, s_rC, F>;
  using ConstSourceFunctionType = ConstDiscreteFunction<V, SGV, s_r, s_rC, F>;

  using RangeSpaceType = SpaceInterface<RGV, r_r, r_rC, F>;
  using RangeFunctionType = DiscreteFunction<V, RGV, r_r, r_rC, F>;
  using ConstRangeFunctionType = ConstDiscreteFunction<V, RGV, r_r, r_rC, F>;

  using ThisType = OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using MatrixOperatorType = MatrixOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using ConstLincombOperatorType = ConstLincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using LincombOperatorType = LincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

  explicit OperatorInterface(const XT::Common::ParameterType& param_type = {})
    : XT::Common::ParametricInterface(param_type)
  {
  }

  virtual ~OperatorInterface() = default;

  /// \name These methods have to be implemented.
  /// \{

  virtual bool linear() const = 0;

  virtual const SourceSpaceType& source_space() const = 0;

  virtual const RangeSpaceType& range_space() const = 0;

  /// \}

  virtual LincombOperatorType operator*(const FieldType& alpha)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, alpha);
    return ret;
  }

  virtual LincombOperatorType operator/(const FieldType& alpha)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1. / alpha);
    return ret;
  }

  virtual ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, 1.);
    return ret;
  }

  virtual ConstLincombOperatorType operator+(const LincombOperatorType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, 1.);
    return ret;
  }

  virtual ConstLincombOperatorType operator+(const ThisType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, 1.);
    return ret;
  }

  virtual LincombOperatorType operator+(LincombOperatorType& other)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, 1.);
    return ret;
  }

  virtual LincombOperatorType operator+(ThisType& other)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, 1.);
    return ret;
  }

  virtual ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, -1.);
    return ret;
  }

  virtual ConstLincombOperatorType operator-(const LincombOperatorType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, -1.);
    return ret;
  }

  virtual ConstLincombOperatorType operator-(const ThisType& other) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, -1.);
    return ret;
  }

  virtual LincombOperatorType operator-(LincombOperatorType& other)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, -1.);
    return ret;
  }

  virtual LincombOperatorType operator-(ThisType& other)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(other, -1.);
    return ret;
  }

  /**
   * Allows the implementation to do preparatory work (i.e., assemble the matrix of a matrix-based operator).
   *
   * \note In general, you have to call this method before calling apply, apply2, apply_inverse or induced_norm!
   */
  virtual ThisType& assemble(const bool /*use_tbb*/ = false)
  {
    return *this;
  }

  /// \name These methods should be implemented and define the functionality of the operator.
  /// \{

  virtual void
  apply(const VectorType& /*source*/, VectorType& /*range*/, const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator cannot be applied!");
  }

  /// \}
  /// \name These methods should be implemented and define the functionality of the operators inverse.
  /// \{

  /**
   * Returns a list of string identifiers for different inversion algorithms, in decending order based on
   * what-we-think-best-ness.
   */
  virtual std::vector<std::string> invert_options() const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator is not invertible!");
    return std::vector<std::string>();
  }

  /**
   * For each inversion algorithm identifier returned by invert_options(), returns an options dict with at least the
   * key "type" set to type
\code
invert_options(some_type).get<std::string>("type") == some_type
\endcode
   * and possibly other key/value pairs.
   */
  virtual XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator is not invertible!");
    return XT::Common::Configuration();
  }

  virtual void apply_inverse(const VectorType& /*range*/,
                             VectorType& /*source*/,
                             const XT::Common::Configuration& /*opts*/,
                             const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator is not invertible!");
  }

  /// \}
  /// \name These methods should be implemented and define the functionality of the operators jacobian.
  /// \{

  /**
   * Same as invert_options(), but concerning the application of the jacobian.
   *
   * \sa invert_options
   */
  virtual std::vector<std::string> jacobian_options() const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator does not provide a jacobian!");
    return std::vector<std::string>();
  }

  /**
   * Same as invert_options(), but concerning the application of the jacobian.
   *
   * \sa invert_options
   */
  virtual XT::Common::Configuration jacobian_options(const std::string& /*type*/) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator does not provide a jacobian!");
    return XT::Common::Configuration();
  }

  /**
   * Either appends suitable functors to the jacobian_op (such that the jacobian of this operator is assembled
   * additively into jacobian_op) or adds the jacobian of this operator to jacobian_op.matrix().
   *
   * \note That you need to call jacobian_op.assemble() to be sure to have the jacobian fully assembled.
   **/
  virtual void jacobian(const VectorType& /*source*/,
                        MatrixOperatorType& /*jacobian_op*/,
                        const XT::Common::Configuration& /*opts*/,
                        const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator does not provide a jacobian!");
  }

  /// \}
  /// \name These apply variants are provided for convenience.
  /// \{

  virtual VectorType apply(const VectorType& source, const XT::Common::Parameter& param = {}) const
  {
    VectorType range(this->range_space().mapper().size(), 0);
    this->apply(source, range, param);
    return range;
  }

  /// \}
  /// \name These apply2 variants are provided for convenience.
  /// \{

  virtual FieldType
  apply2(const VectorType& range, const VectorType& source, const XT::Common::Parameter& param = {}) const
  {
    return range.dot(this->apply(source, param));
  }

  /// \}
  /// \name These apply_invers variants are provided for convenience.
  /// \{

  virtual void apply_inverse(const VectorType& range,
                             VectorType& source,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options(type), param);
  }

  virtual void apply_inverse(const VectorType& range, VectorType& source, const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options().at(0), param);
  }

  virtual VectorType apply_inverse(const VectorType& range,
                                   const XT::Common::Configuration& opts,
                                   const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, opts, param);
    return source;
  }

  virtual VectorType
  apply_inverse(const VectorType& range, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, type, param);
    return source;
  }

  virtual VectorType apply_inverse(const VectorType& range, const XT::Common::Parameter& param = {}) const
  {
    VectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, param);
    return source;
  }

  /// \}
  /// \name These jacobian variants are provided for convenience.
  /// \{

  virtual void jacobian(const VectorType& source,
                        MatrixOperatorType& jacobian_op,
                        const std::string& type,
                        const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, jacobian_op, this->jacobian_options(type), param);
  }

  virtual void
  jacobian(const VectorType& source, MatrixOperatorType& jacobian_op, const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, jacobian_op, this->jacobian_options().at(0), param);
  }

  virtual MatrixOperatorType jacobian(const VectorType& source,
                                      const XT::Common::Configuration& opts,
                                      const XT::Common::Parameter& param = {}) const
  {
    MatrixOperatorType jacobian_op(this->source_space().grid_view(),
                                   this->source_space(),
                                   this->range_space(),
                                   make_element_and_intersection_sparsity_pattern(
                                       this->range_space(), this->source_space(), this->source_space().grid_view()));
    this->jacobian(source, jacobian_op, opts, param);
    return jacobian_op;
  } // ... jacobian(...)

  virtual MatrixOperatorType
  jacobian(const VectorType& source, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    MatrixOperatorType jacobian_op(this->source_space().grid_view(),
                                   this->source_space(),
                                   this->range_space(),
                                   make_element_and_intersection_sparsity_pattern(
                                       this->range_space(), this->source_space(), this->source_space().grid_view()));
    this->jacobian(source, jacobian_op, type, param);
    return jacobian_op;
  } // ... jacobian(...)

  virtual MatrixOperatorType jacobian(const VectorType& source, const XT::Common::Parameter& param = {}) const
  {
    MatrixOperatorType jacobian_op(this->source_space().grid_view(),
                                   this->source_space(),
                                   this->range_space(),
                                   make_element_and_intersection_sparsity_pattern(
                                       this->range_space(), this->source_space(), this->source_space().grid_view()));
    this->jacobian(source, jacobian_op, param);
    return jacobian_op;
  } // ... jacobian(...)

  /// \}
  /// \name These induced_norm variants are provided for convenience.
  /// \{

  template <class ParameterType_,
            typename = /* Only enable this method, if */
            typename std::enable_if</* param is the same as XT::Common::Parameter */ (
                                        std::is_same<ParameterType_, XT::Common::Parameter>::value)
                                    && /* and the vector spaces defined by SourceSpaceType/VectorType and */
                                    /* RangeSpaceType/VectorType coincide. */ (
                                        std::is_same<V, V>::value&& std::is_same<SGV, RGV>::value && (s_r == r_r)
                                        && (s_rC == r_rC))>::type>
  FieldType induced_norm(const VectorType& range, const ParameterType_& param = {}) const
  {
    using std::sqrt;
    return sqrt(this->apply2(range, range, param));
  }

  /// \}
  /// \name For each method above which accepts vectors we provide a similar variant for convenience which accepts
  ///       functions (which simply extracts the vectors and calls the respectiv method).
  /// \{

  virtual void
  apply(const ConstSourceFunctionType& source, RangeFunctionType& range, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    this->apply(source.dofs().vector(), range.dofs().vector(), param);
  }

  virtual RangeFunctionType apply(const ConstSourceFunctionType& source, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return RangeFunctionType(this->range_space(), this->apply(source.dofs().vector(), param));
  }

  virtual FieldType apply2(const SourceFunctionType& range,
                           const RangeFunctionType& source,
                           const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    return this->apply2(range.dofs().vector(), source.dofs().vector(), param);
  }

  virtual void apply_inverse(const RangeFunctionType& range,
                             SourceFunctionType& source,
                             const XT::Common::Configuration& opts,
                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    this->apply_inverse(range.dofs().vector(), source.dofs().vector(), opts, param);
  }

  virtual void apply_inverse(const RangeFunctionType& range,
                             SourceFunctionType& source,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    this->apply_inverse(range.dofs().vector(), source.dofs().vector(), type, param);
  }

  virtual void apply_inverse(const RangeFunctionType& range,
                             SourceFunctionType& source,
                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    this->apply_inverse(range.dofs().vector(), source.dofs().vector(), param);
  }

  virtual SourceFunctionType apply_inverse(const RangeFunctionType& range,
                                           const XT::Common::Configuration& opts,
                                           const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    return SourceFunctionType(this->source_space(), this->apply_inverse(range.dofs().vector(), opts, param));
  }

  virtual SourceFunctionType
  apply_inverse(const RangeFunctionType& range, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    return SourceFunctionType(this->source_space(), this->apply_inverse(range.dofs().vector(), type, param));
  }

  virtual SourceFunctionType apply_inverse(const RangeFunctionType& range,
                                           const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    return SourceFunctionType(this->source_space(), this->apply_inverse(range.dofs().vector(), param));
  }

  virtual void jacobian(const SourceFunctionType& source,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Configuration& opts,
                        const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), jacobian_op, opts, param);
  }

  virtual void jacobian(const SourceFunctionType& source,
                        MatrixOperatorType& jacobian_op,
                        const std::string& type,
                        const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), jacobian_op, type, param);
  }

  virtual void jacobian(const SourceFunctionType& source,
                        MatrixOperatorType& jacobian_op,
                        const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), jacobian_op, param);
  }

  virtual MatrixOperatorType jacobian(const SourceFunctionType& source,
                                      const XT::Common::Configuration& opts,
                                      const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), opts, param);
  }

  virtual MatrixOperatorType
  jacobian(const SourceFunctionType& source, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), type, param);
  }

  virtual MatrixOperatorType jacobian(const SourceFunctionType& source, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), param);
  }

  template <class ParameterType_,
            typename = /* Only enable this method, if */
            typename std::enable_if</* param is the same as XT::Common::Parameter */ (
                                        std::is_same<ParameterType_, XT::Common::Parameter>::value)
                                    && /* and the vector spaces defined by SourceSpaceType/VectorType and */
                                    /* RangeSpaceType/VectorType coincide. */ (
                                        std::is_same<V, V>::value&& std::is_same<SGV, RGV>::value && (s_r == r_r)
                                        && (s_rC == r_rC))>::type>
  FieldType induced_norm(const RangeFunctionType& range, const ParameterType_& param = {}) const
  {
    DUNE_THROW_IF(!this->range_space().contains(range),
                  Exceptions::operator_error,
                  "this->range_space() = " << this->range_space() << "\n   range.space() = " << range.space());
    return this->induced_norm(range.dofs().vector(), param);
  }

  /// \}
}; // class OperatorInterface


} // namespace GDT
} // namespace Dune

#include "lincomb.hh"
#include "matrix-based.hh"

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
