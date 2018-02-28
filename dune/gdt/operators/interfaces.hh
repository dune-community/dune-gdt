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
 * Consider (discrete) spaces V_h and W_h, and a field K, this interface models
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
 * The functions v_h \in V_h (and w_h \in W_h), to which the operator can be applied to, are modelled by
 * ConstDiscreteFunction or DiscreteFunciton with an appropriate vector type derived from XT::LA::VectorInterface
 * modelled by SourceVector (RangeVector).
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

  using RangeSpaceType = SpaceInterface<RGV, r_r, r_rC, RF>;
  using RangeVectorType = RangeVector;

  using FieldType = Field;

  using ThisType = OperatorInterface<SV, SGV, s_r, s_rC, SF, F, r_r, r_rC, RF, RGV, RV>;

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

  /// \name These methods should be implemented.
  /// \{

  /**
   * \note Should raise an operator_error if this operator is only a two-form.
   **/
  virtual void apply(const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& /*source*/,
                     DiscreteFunction<RV, RGV, r_r, r_rC, RF>& /*range*/,
                     const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error,
               "Either this two-form is not an operator or the implementor forgot to override apply!");
  }

  /**
   * Returns a list of string identifiers for different inversion algorithms, in decending order based on
   * what-we-think-best-ness.
   *
   * \note Should raise an operator_error if this operator is not invertible, as opposed to an empty list.
   */
  virtual std::vector<std::string> invert_options() const
  {
    DUNE_THROW(Exceptions::operator_error,
               "Either this operator is not invertible or the implementor forgot to override invert_options!");
  }

  /**
   * For each inversion algorithm identifier returned by invert_options(), returns an options dict with at least the
   * key "type" set to type
\code
invert_options(some_type).get<std::string>("type") == some_type
\endcode
   * and possible other key/value pairs.
   */
  virtual XT::Common::Configuration invert_options(const std::string& /*type*/) const
  {
    DUNE_THROW(Exceptions::operator_error,
               "Either this operator is not invertible or the implementor forgot to override invert_options!");
  }

  /**
   * Should raise an operator_error if this operator is not invertible.
   **/
  virtual void apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& /*range*/,
                             DiscreteFunction<SV, SGV, s_r, s_rC, SF>& /*source*/,
                             const XT::Common::Configuration& /*opts*/,
                             const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error,
               "Either this operator is not invertible or the implementor forgot to override apply_inverse!");
  }

  /**
   * Same as invert_options(), but concerning the application of the jacobian.
   *
   * \sa invert_options
   */
  virtual std::vector<std::string> jacobian_options() const
  {
    DUNE_THROW(
        Exceptions::operator_error,
        "Either this operator does not provide a jacobian or the implementor forgot to override jacobian_options!");
  }

  /**
   * Same as invert_options(), but concerning the application of the jacobian.
   *
   * \sa invert_options
   */
  virtual XT::Common::Configuration jacobian_options(const std::string& /*type*/) const
  {
    DUNE_THROW(
        Exceptions::operator_error,
        "Either this operator does not provide a jacobian or the implementor forgot to override jacobian_options!");
  }

  /**
   * Should raise an operator_error if this operator does not provide a jacobian.
   **/
  virtual std::shared_ptr<ThisType> jacobian(const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& /*source*/,
                                             const XT::Common::Configuration& /*opts*/,
                                             const XT::Common::Parameter& /*param*/ = {}) const
  {
    DUNE_THROW(Exceptions::operator_error,
               "Either this operator does not provide a jacobian or the implementor forgot to override jacobian!");
  }

  /// \}
  /// \name These methods are provided for convenience.
  /// \{

  virtual DiscreteFunction<RV, RGV, r_r, r_rC, RF> apply(const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                                                         const XT::Common::Parameter& param = {}) const
  {
    DiscreteFunction<RV, RGV, r_r, r_rC, RF> range(this->range_space());
    this->apply(source, range, param);
    return range;
  }

  virtual FieldType apply2(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                           const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                           const XT::Common::Parameter& param = {}) const
  {
    return range.dofs().vector().dot(this->apply(source, param).dofs().vector());
  }

  virtual void apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                             DiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                             const std::string& type,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options(type), param);
  }

  virtual void apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                             DiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                             const XT::Common::Parameter& param = {}) const
  {
    this->apply_inverse(range, source, this->invert_options().at(0), param);
  }

  virtual DiscreteFunction<SV, SGV, s_r, s_rC, SF>
  apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const
  {
    DiscreteFunction<SV, SGV, s_r, s_rC, SF> source(this->source_space());
    this->apply_inverse(range, source, opts, param);
    return source;
  }

  virtual DiscreteFunction<SV, SGV, s_r, s_rC, SF>
  apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                const std::string& type,
                const XT::Common::Parameter& param = {}) const
  {
    DiscreteFunction<SV, SGV, s_r, s_rC, SF> source(this->source_space());
    this->apply_inverse(range, source, type, param);
    return source;
  }

  virtual DiscreteFunction<SV, SGV, s_r, s_rC, SF>
  apply_inverse(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                const XT::Common::Parameter& param = {}) const
  {
    DiscreteFunction<SV, SGV, s_r, s_rC, SF> source(this->source_space());
    this->apply_inverse(range, source, param);
    return source;
  }

  virtual std::shared_ptr<ThisType> jacobian(const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                                             const std::string& type,
                                             const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, this->jacobian_options(type), param);
  }

  virtual std::shared_ptr<ThisType> jacobian(const ConstDiscreteFunction<SV, SGV, s_r, s_rC, SF>& source,
                                             const XT::Common::Parameter& param = {}) const
  {
    return this->jacobian(source, this->jacobian_options().at(0), param);
  }

  template <class RV_,
            class RGV_,
            size_t r_r_,
            size_t r_rC_,
            class RF_,
            typename = /* Only enable this method, if */
            typename std::enable_if</* range is in the vector space defined by RangeSpaceType/RangeVectorType */ (
                                        std::is_same<RV_, RV>::value&& std::is_same<RGV_, RGV>::value && (r_r_ == r_r)
                                        && (r_rC_ == r_rC)
                                        && std::is_same<RF_, RF>::value)
                                    && /* and the vector spaces defined by SourceSpaceType/SourceVectorType and */
                                    /* RangeSpaceType/RangeVectorType coincide. */ (
                                        std::is_same<SV, RV>::value&& std::is_same<SGV, RGV>::value && (s_r == r_r)
                                        && (s_rC == r_rC)
                                        && std::is_same<SF, RF>::value)>::type>
  FieldType induced_norm(const ConstDiscreteFunction<RV, RGV, r_r, r_rC, RF>& range,
                         const XT::Common::Parameter& param = {}) const
  {
    return std::sqrt(this->apply2(range, range, param));
  }

  /// \}
}; // class OperatorInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
