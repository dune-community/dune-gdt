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
#include <dune/xt/common/timedlogging.hh>
#include <dune/xt/common/type_traits.hh>
#include <dune/xt/la/exceptions.hh>
#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/la/solver.hh>
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
class ConstantOperator;

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
class OperatorInterface
  : internal::AssertArgumentsOfOperatorInterface<Matrix,
                                                 SourceGridView,
                                                 source_dim,
                                                 source_dim_cols,
                                                 range_dim,
                                                 range_dim_cols,
                                                 RangeGridView>
  , public XT::Common::ParametricInterface
{
public:
  using MatrixType = Matrix;
  using VectorType = XT::LA::vector_t<MatrixType>;
  using FieldType = typename MatrixType::ScalarType;

  using M = Matrix;
  using V = VectorType;
  using F = FieldType;

  using SGV = SourceGridView;
  using E = XT::Grid::extract_entity_t<SGV>;
  static const constexpr size_t s_r = source_dim;
  static const constexpr size_t s_rC = source_dim_cols;

  using RGV = RangeGridView;
  static const constexpr size_t r_r = range_dim;
  static const constexpr size_t r_rC = range_dim_cols;

  using SourceSpaceType = SpaceInterface<SGV, s_r, s_rC, F>;
  using SourceFunctionInterfaceType = XT::Functions::GridFunctionInterface<E, s_r, s_rC, F>;
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
  {}

  virtual ~OperatorInterface() = default;

  /// \name These methods have to be implemented.
  /// \{

  virtual bool linear() const = 0;

  virtual const SourceSpaceType& source_space() const = 0;

  virtual const RangeSpaceType& range_space() const = 0;

  /// \}

  virtual ConstLincombOperatorType operator*(const FieldType& alpha) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, alpha);
    return ret;
  }

  virtual LincombOperatorType operator*(const FieldType& alpha)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, alpha);
    return ret;
  }

  virtual ConstLincombOperatorType operator/(const FieldType& alpha) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1. / alpha);
    return ret;
  }

  virtual LincombOperatorType operator/(const FieldType& alpha)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1. / alpha);
    return ret;
  }

  /// \name const operator+ variants
  /// \{

  virtual ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const
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

  /**
   * \note vector is interpreted as a ConstantOperator
   * \sa ConstantOperator
   */
  virtual ConstLincombOperatorType operator+(const VectorType& vector) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            1.);
    return ret;
  }

  /// \}
  /// \name mutable operator+ variants
  /// \{

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

  /**
   * \note vector is interpreted as a ConstantOperator
   * \sa ConstantOperator
   */
  virtual LincombOperatorType operator+(const VectorType& vector)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            1.);
    return ret;
  }

  /// \}
  /// \name const operator- variants
  /// \{

  virtual ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const
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

  /**
   * \note vector is interpreted as a ConstantOperator
   * \sa ConstantOperator
   */
  virtual ConstLincombOperatorType operator-(const VectorType& vector) const
  {
    ConstLincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            -1.);
    return ret;
  }

  /// \}
  /// \name mutable operator- variants between arbitrary operators
  /// \{

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
   * \note vector is interpreted as a ConstantOperator
   * \sa ConstantOperator
   */
  virtual LincombOperatorType operator-(const VectorType& vector)
  {
    LincombOperatorType ret(this->source_space(), this->range_space());
    ret.add(*this, 1.);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            -1.);
    return ret;
  }

  /// \}

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
    return {"newton"};
  }

  /**
   * For each inversion algorithm identifier returned by invert_options(), returns an options dict with at least the
   * key "type" set to type
\code
invert_options(some_type).get<std::string>("type") == some_type
\endcode
   * and possibly other key/value pairs.
   *
   * \todo Allow to pass jacobian options as subcfg in newton.
   */
  virtual XT::Common::Configuration invert_options(const std::string& type) const
  {
    if (type == "newton") {
      return {{"type", type}, {"precision", "1e-7"}, {"max_iter", "100"}, {"max_dampening_iter", "1000"}};
    } else
      DUNE_THROW(Exceptions::operator_error, "type = " << type);
  } // ... invert_options(...)

  /**
   * Given source is used as initial guess.
   *
   * \todo Allow to pass jacobian options as subcfg in newton.
   **/
  virtual void apply_inverse(const VectorType& range,
                             VectorType& source,
                             const XT::Common::Configuration& opts,
                             const XT::Common::Parameter& param = {}) const
  {
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, "opts = " << opts);
    const auto type = opts.get<std::string>("type");
    const XT::Common::Configuration default_opts = this->invert_options(type);
    auto logger = XT::Common::TimedLogger().get("gdt.operator.apply_inverse");
    if (type == "newton") {
      // some preparations
      auto residual_op = *this - range;
      auto residual = range.copy();
      auto update = source.copy();
      auto candidate = source.copy();
      // one matrix for all jacobians
      MatrixOperatorType jacobian_op(this->source_space().grid_view(),
                                     this->source_space(),
                                     this->range_space(),
                                     make_element_and_intersection_sparsity_pattern(
                                         this->source_space(), this->range_space(), this->source_space().grid_view()));
      XT::LA::Solver<M> jacobian_solver(jacobian_op.matrix());
      const auto precision = opts.get("precision", default_opts.get<double>("precision"));
      const auto max_iter = opts.get("max_iter", default_opts.get<size_t>("max_iter"));
      const auto max_dampening_iter = opts.get("max_dampening_iter", default_opts.get<size_t>("max_iter"));
      size_t l = 0;
      Timer timer;
      while (true) {
        timer.reset();
        logger.debug() << "l = " << l << ": computing residual ... " << std::flush;
        residual_op.apply(source, residual, param);
        auto res = residual.l2_norm();
        logger.debug() << "took " << timer.elapsed() << "s, |residual|_l2 = " << res << std::endl;
        if (res < precision) {
          logger.debug() << "       residual below tolerance, succeeded!" << std::endl;
          break;
        }
        DUNE_THROW_IF(l >= max_iter,
                      Exceptions::operator_error,
                      "max iterations reached!\n|residual|_l2 = " << res << "\nopts:\n"
                                                                  << opts);
        logger.debug() << "       computing jacobi matrix ... " << std::flush;
        timer.reset();
        jacobian_op.matrix() *= 0.;
        residual_op.jacobian(source, jacobian_op, {{"type", residual_op.jacobian_options().at(0)}}, param);
        jacobian_op.walk(/*use_tbb=*/true);
        logger.debug() << "took " << timer.elapsed() << "s"
                       << "\n"
                       << "       solving for defect ... " << std::flush;
        timer.reset();
        residual *= -1.;
        update = source; // <- initial guess for the linear solver
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
        logger.debug() << "took " << timer.elapsed() << "s";
        if (tried_linear_solvers.size() > 1)
          logger.debug() << ", took " << tried_linear_solvers.size() << " attempts with different linear solvers";
        logger.debug() << "\n       computing update ... " << std::flush;
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
          candidate = source + update * lambda;
          residual_op.apply(candidate, residual, param);
          candidate_res = residual.l2_norm();
          lambda /= 2;
          k += 1;
        }
        source = candidate;
        logger.debug() << "took " << timer.elapsed() << "s and a dampening of " << 2 * lambda << std::endl;
        l += 1;
      }
    } else
      DUNE_THROW(Exceptions::operator_error, "type = " << type);
  } // ... apply_inverse(...)

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
    VectorType range(this->range_space().mapper().size(), 0.);
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

  template <
      class ParameterType_,
      typename = /* Only enable this method, if */
      typename std::enable_if<
          /* param is the same as XT::Common::Parameter */ (std::is_same<ParameterType_, XT::Common::Parameter>::value)
          && /* and the vector spaces defined by SourceSpaceType/VectorType and */
          /* RangeSpaceType/VectorType coincide. */ (std::is_same<V, V>::value&& std::is_same<SGV, RGV>::value
                                                     && (s_r == r_r) && (s_rC == r_rC))>::type>
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

  virtual FieldType apply2(const RangeFunctionType& range,
                           const SourceFunctionType& source,
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

  template <
      class ParameterType_,
      typename = /* Only enable this method, if */
      typename std::enable_if<
          /* param is the same as XT::Common::Parameter */ (std::is_same<ParameterType_, XT::Common::Parameter>::value)
          && /* and the vector spaces defined by SourceSpaceType/VectorType and */
          /* RangeSpaceType/VectorType coincide. */ (std::is_same<V, V>::value&& std::is_same<SGV, RGV>::value
                                                     && (s_r == r_r) && (s_rC == r_rC))>::type>
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

#include "constant.hh"
#include "lincomb.hh"
#include "matrix-based.hh"

#endif // DUNE_GDT_OPERATORS_INTERFACES_HH
