// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
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

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Interface for operators (and two-forms).
 *
 * Considering (discrete) spaces V_h and W_h, and a field K, this interface models
 *
 * - operators A: V_h -> W_h and
 *
 * - two-forms B: W_h x V_h -> K.
 *
 * The source of the operator V_h is the (discrete) space
 *
 *   V_h := \{ v_h: \tau_h^s -> SF^{s_r \times s_rC} \}
 *
 * (modelled by SpaceInterface) of functions mapping from a partition \tau_h^s (modelled by SourceGridView) of a
 * physical domain to a (possibly vector- or matrix-valued) vector space SF^{s_r \times s_rC} (modelled by SourceField,
 * source_dim and source_dim_cols). The range of the operator W_h (identified with its dual, since we are in the
 * discrete setting) is the (discrete) space
 *
 *   W_h := \{ w_h: \tau_h^r -> RF^{r_r \times r_rC} \}
 *
 * (modelled by SpaceInterface) of functions mapping from a partition \tau_h^r (modelled by RangeGridView) of a
 * physical domain to a (possibly vector- or matrix-valued) vector space RF^{r_r \times r_rC} (modelled by RangeField,
 * range_dim and range_dim_cols).
 *
 * The functions v_h \in V_h (and w_h \in W_h), to which the operator can be applied to, are represented by their DoF
 * vectors (an appropriate vector type derived from XT::LA::VectorInterface modelled by SourceVector (RangeVector)),
 * which are then interpreted as discrete functions in the respective source_space or range_space.
 *
 * The field K of the interpretation of the operator as a two-form (see for instance the default implementation of
 * apply2()) is modelled by Field.
 */
template <class SourceVector,
          class SourceGridView,
          size_t source_dim = 1,
          size_t source_dim_cols = 1,
          class SourceField = double,
          class Field = double,
          size_t range_dim = source_dim,
          size_t range_dim_cols = source_dim_cols,
          class RangeField = double,
          class RangeGridView = SourceGridView,
          class RangeVector = SourceVector>
class OperatorInterface : public XT::Common::ParametricInterface
{
  static_assert(XT::LA::is_vector<SourceVector>::value, "");
  static_assert(XT::LA::is_vector<RangeVector>::value, "");

public:
  using SV = SourceVector;
  using SGV = SourceGridView;
  static const constexpr size_t s_r = source_dim;
  static const constexpr size_t s_rC = source_dim_cols;
  using SF = SourceField;

  using RV = RangeVector;
  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_dim;
  static const constexpr size_t r_rC = range_dim_cols;
  using RF = RangeField;

  using F = Field;

  using SourceSpaceType = SpaceInterface<SGV, s_r, s_rC, SF>;
  using SourceVectorType = SourceVector;
  using SourceFunctionType = DiscreteFunction<SourceVector, SGV, s_r, s_rC, SF>;
  using ConstSourceFunctionType = ConstDiscreteFunction<SourceVector, SGV, s_r, s_rC, SF>;

  using RangeSpaceType = SpaceInterface<RGV, r_r, r_rC, RF>;
  using RangeVectorType = RangeVector;
  using RangeFunctionType = DiscreteFunction<RangeVector, RGV, r_r, r_rC, RF>;
  using ConstRangeFunctionType = ConstDiscreteFunction<RangeVector, RGV, r_r, r_rC, RF>;

  using FieldType = Field;

  using ThisType = OperatorInterface<SV, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, RV>;

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

  /**
   * Allows the implementation to do preparatory work (i.e., assemble the matrix of a matrix-based linear operator).
   *
   * \note In general, you have to call this method before calling apply, apply2, apply_inverse or induced_norm!
   */
  virtual void assemble(const bool /*use_tbb*/ = 0)
  {
  }

  /// \name These methods should be implemented and define the functionality of the operator.
  /// \{

  virtual void apply(const SourceVectorType& /*source*/,
                     RangeVectorType& /*range*/,
                     const XT::Common::Parameter& /*param*/ = {}) const
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

  virtual void apply_inverse(const RangeVectorType& /*range*/,
                             SourceVectorType& /*source*/,
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

  virtual std::shared_ptr<ThisType> jacobian(const SourceVectorType& /*source*/,
                                             const XT::Common::Configuration& /*opts*/,
                                             const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error, "This operator does not provide a jacobian!");
    return nullptr;
  }

  /// \}
  /// \name These apply variants are provided for convenience.
  /// \{

  virtual RangeVectorType apply(const SourceVectorType& source, const XT::Common::Parameter& param = {}) const
  {
    RangeVectorType range(this->range_space().mapper().size(), 0);
    this->apply(source, range, param);
    return range;
  }

  /// \}
  /// \name These apply2 variants are provided for convenience.
  /// \{

  virtual FieldType
  apply2(const SourceVectorType& range, const RangeVectorType& source, const XT::Common::Parameter& param = {}) const
  {
    return range.dot(this->apply(source, param));
  }

  /// \}
  /// \name These apply_invers variants are provided for convenience.
  /// \{

  virtual void apply_inverse(const RangeVectorType& range,
                             SourceVectorType& source,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options(type), param);
  }

  virtual void
  apply_inverse(const RangeVectorType& range, SourceVectorType& source, const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options().at(0), param);
  }

  virtual SourceVectorType apply_inverse(const RangeVectorType& range,
                                         const XT::Common::Configuration& opts,
                                         const XT::Common::Parameter& param = {}) const
  {
    SourceVectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, opts, param);
    return source;
  }

  virtual SourceVectorType
  apply_inverse(const RangeVectorType& range, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    SourceVectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, type, param);
    return source;
  }

  virtual SourceVectorType apply_inverse(const RangeVectorType& range, const XT::Common::Parameter& param = {}) const
  {
    SourceVectorType source(this->source_space().mapper().size());
    this->apply_inverse(range, source, param);
    return source;
  }

  /// \}
  /// \name These jacobian variants are provided for convenience.
  /// \{

  virtual std::shared_ptr<ThisType>
  jacobian(const SourceVectorType& source, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, this->jacobian_options(type), param);
  }

  virtual std::shared_ptr<ThisType> jacobian(const SourceVectorType& source,
                                             const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, this->jacobian_options().at(0), param);
  }

  /// \}
  /// \name These induced_norm variants are provided for convenience.
  /// \{

  template <class ParameterType_,
            typename = /* Only enable this method, if */
            typename std::enable_if</* param is the same as XT::Common::Parameter */ (
                                        std::is_same<ParameterType_, XT::Common::Parameter>::value)
                                    && /* and the vector spaces defined by SourceSpaceType/SourceVectorType and */
                                    /* RangeSpaceType/RangeVectorType coincide. */ (
                                        std::is_same<SV, RV>::value&& std::is_same<SGV, RGV>::value && (s_r == r_r)
                                        && (s_rC == r_rC)
                                        && std::is_same<SF, RF>::value)>::type>
  FieldType induced_norm(const RangeVectorType& range, const ParameterType_& param = {}) const
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

  virtual std::shared_ptr<ThisType> jacobian(const SourceFunctionType& source,
                                             const XT::Common::Configuration& opts,
                                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), opts, param);
  }

  virtual std::shared_ptr<ThisType>
  jacobian(const SourceFunctionType& source, const std::string& type, const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!this->source_space().contains(source),
                  Exceptions::operator_error,
                  "this->source_space() = " << this->source_space() << "\n   source.space() = " << source.space());
    return this->jacobian(source.dofs().vector(), type, param);
  }

  virtual std::shared_ptr<ThisType> jacobian(const SourceFunctionType& source,
                                             const XT::Common::Parameter& param = {}) const
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
                                    && /* and the vector spaces defined by SourceSpaceType/SourceVectorType and */
                                    /* RangeSpaceType/RangeVectorType coincide. */ (
                                        std::is_same<SV, RV>::value&& std::is_same<SGV, RGV>::value && (s_r == r_r)
                                        && (s_rC == r_rC)
                                        && std::is_same<SF, RF>::value)>::type>
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

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
