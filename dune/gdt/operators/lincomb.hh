// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_LINCOMB_HH
#define DUNE_GDT_OPERATORS_LINCOMB_HH

#include <vector>

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>

#include "interfaces.hh"

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
class ConstLincombOperator : public OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
{
public:
  using ThisType = ConstLincombOperator;
  using BaseType = OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::ConstLincombOperatorType;
  using typename BaseType::FieldType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using OperatorType = BaseType;

  ConstLincombOperator(const AssemblyGridViewType& assembly_grid_vw,
                       const SourceSpaceType& src_space,
                       const RangeSpaceType& rng_space,
                       const std::string& logging_prefix = "",
                       const std::array<bool, 3>& logging_enabled = XT::Common::default_logger_state())
    : BaseType({}, logging_prefix.empty() ? "ConstLincombOperator" : logging_prefix, logging_enabled)
    , assembly_grid_view_(assembly_grid_vw)
    , source_space_(src_space)
    , range_space_(rng_space)
  {
    LOG_(debug) << "ConstLincombOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &src_space
                << ", range_space=" << &rng_space << ")" << std::endl;
  }

  ConstLincombOperator(const ThisType& other) = default;

  ConstLincombOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;
  using BaseType::jacobian;
  using BaseType::operator+;
  using BaseType::operator-;
  using BaseType::operator*;
  using BaseType::operator/;

  /// \name Required by OperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const final
  {
    return range_space_;
  }

  bool linear() const final
  {
    for (const auto& op : const_ops_)
      if (!op.access().linear())
        return false;
    return true;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const final
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    LOG_(info) << "applying " << this->num_ops() << " operators ..." << std::endl;
    range_vector.set_all(0);
    auto tmp_range_vector = range_vector;
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      this->op(ii).apply(source_function, tmp_range_vector, param);
      tmp_range_vector *= this->coeff(ii);
      range_vector += tmp_range_vector;
    }
  } // ... apply(...)

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const final
  {
    return source_space_;
  }

  const AssemblyGridViewType& assembly_grid_view() const final
  {
    return assembly_grid_view_;
  }

  void
  apply(const VectorType& source_vector, VectorType& range_vector, const XT::Common::Parameter& param = {}) const final
  {
    LOG_(debug) << "apply(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    LOG_(info) << "applying " << this->num_ops() << " operators ..." << std::endl;
    range_vector.set_all(0);
    auto tmp_range_vector = range_vector;
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      this->op(ii).apply(source_vector, tmp_range_vector, param);
      tmp_range_vector *= this->coeff(ii);
      range_vector += tmp_range_vector;
    }
  } // ... apply(...)

protected:
  std::vector<XT::Common::Configuration> all_jacobian_options() const final
  {
    std::vector<XT::Common::Configuration> ret(1);
    auto& cfg = ret[0];
    cfg["type"] = "lincomb";
    using XT::Common::to_string;
    for (size_t ii = 0; ii < this->num_ops(); ++ii)
      cfg.add(this->op(ii).jacobian_options(this->op(ii).jacobian_options().at(0)), "op." + to_string(ii));
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      const auto back_ii = ssize_t(this->num_ops()) - 1 - ssize_t(ii);
      cfg.add(this->op(back_ii).jacobian_options(this->op(back_ii).jacobian_options().at(0)),
              "back_op." + to_string(back_ii));
    }
    return ret;
  } // ... all_jacobian_options(...)

public:
  void jacobian(const VectorType& source_vector,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override
  {

    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm()
                << ",\n   opts=" << print(opts, {{"oneline", "true"}}) << ",\n   param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_source(source_vector);
    this->assert_jacobian_opts(opts); // ensures that "lincomb" is requested
    using XT::Common::to_string;
    LOG_(info) << "appending " << this->num_ops() << " jacobians ..." << std::endl;
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      // parse opts
      XT::Common::Configuration op_opts;
      if (opts.has_sub("op." + to_string(ii) + ".opts"))
        op_opts = opts.sub("op." + to_string(ii));
      const auto back_ii = ssize_t(this->num_ops()) - 1 - ssize_t(ii);
      if (opts.has_sub("back_op." + to_string(back_ii) + ".opts")) {
        auto back_op_opts = opts.sub("back_op." + to_string(back_ii));
        DUNE_THROW_IF(!opts.empty() && !back_op_opts.empty() && (op_opts != back_op_opts),
                      Exceptions::operator_error,
                      "Conflicting opts specified for operator " << ii << " (once as op." << ii << ", onse as back_op."
                                                                 << back_ii << ", see below)!"
                                                                 << "\n\n"
                                                                 << opts);
      }
      if (op_opts.empty())
        op_opts["type"] = this->op(ii).jacobian_options().at(0);
      // save curent scaling
      const auto scaling = jacobian_op.scaling;
      // add op
      jacobian_op.scaling *= this->coeff(ii);
      LOG_(debug) << "   backing up scaling = " << scaling << "\n   this->coeff(ii) = " << this->coeff(ii)
                  << "\n   jacobian_op.scaling = " << jacobian_op.scaling << std::endl;
      this->op(ii).jacobian(source_vector, jacobian_op, op_opts, param);
      // restore scaling
      LOG_(debug) << "   restoring jacobian_op.scaling = " << scaling << std::endl;
      jacobian_op.scaling = scaling;
    }
  } // ... jacobian(...)

  /// \}
  /// \name These methods allow access to the summands of the linear combination
  /// \{

  size_t num_ops() const
  {
    return const_ops_.size();
  }

  const OperatorType& op(const size_t ii) const
  {
    DUNE_THROW_IF(ii >= this->num_ops(),
                  Exceptions::operator_error,
                  "ii = " << ii << "\n   this->num_ops() = " << this->num_ops());
    return const_ops_[ii].access();
  }

  const FieldType& coeff(const size_t ii) const
  {
    DUNE_THROW_IF(ii >= this->num_ops(),
                  Exceptions::operator_error,
                  "ii = " << ii << "\n   this->num_ops() = " << this->num_ops());
    return coeffs_[ii];
  }

  /// \}
  /// \name These methods allow to add a summand to the linear combination
  /// \{

  void add(const OperatorType& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op.logger.state);
    LOG_(debug) << "add(const_op_ref=" << &op << ", coeff=" << coeff << ")" << std::endl;
    this->extend_parameter_type(op.parameter_type());
    const_ops_.emplace_back(op);
    coeffs_.emplace_back(coeff);
  }

  void add(OperatorType*&& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op->logger.state);
    LOG_(debug) << "add(op_ptr=" << op << ", coeff=" << coeff << ")" << std::endl;
    this->extend_parameter_type(op->parameter_type());
    keep_alive_.emplace_back(std::move(op));
    const_ops_.emplace_back(*keep_alive_.back());
    coeffs_.emplace_back(coeff);
  }

  void add(const ThisType& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op.logger.state);
    LOG_(debug) << "add(const_lincomb_op_ref=" << &op << ", coeff=" << coeff << ")" << std::endl;
    this->extend_parameter_type(op.parameter_type());
    // only adding op itself would lead to segfaults if op is a temporary, so we need to
    // - keep all those operators alive that op cared about
    for (auto shrd_ptr : op.keep_alive_)
      this->keep_alive_.emplace_back(shrd_ptr);
    // - take over ops linear decomposition
    for (size_t ii = 0; ii < op.num_ops(); ++ii) {
      LOG_(debug) << "  adding op=" << &(op.const_ops_[ii]) << ", coeff=" << coeff << std::endl;
      const_ops_.emplace_back(op.const_ops_[ii]);
      coeffs_.emplace_back(coeff * op.coeffs_[ii]);
    }
  } // ... add(...)

  // we need this, otherwise add(OperatorType*&&) would be used
  void add(ThisType*&& op, const FieldType& coeff = 1.)
  {
    this->add(*op, coeff);
  }

  /// \}
  /// \name These numeric operators extend the ones from OperatorInterface
  /// \{

  ThisType& operator*=(const FieldType& alpha)
  {
    LOG_(debug) << "operator*=(alpha=" << alpha << ")" << std::endl;
    for (auto& coeff : coeffs_)
      coeff *= alpha;
    return *this;
  }

  ThisType& operator/=(const FieldType& alpha)
  {
    LOG_(debug) << "operator/=(alpha=" << alpha << ")" << std::endl;
    for (auto& coeff : coeffs_)
      coeff /= alpha;
    return *this;
  }

  ThisType& operator+=(const BaseType& other)
  {
    LOG_(debug) << "operator+=(other_op=" << &other << ")" << std::endl;
    this->add(other);
    return *this;
  }

  ThisType& operator+=(const ThisType& other)
  {
    this->logger.state_or(other.logger.state);
    LOG_(debug) << "operator+=(other_const_lincomb_op=" << &other << ")" << std::endl;
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      LOG_(debug) << "  adding op=" << &(other.const_ops_[ii]) << ", coeff=" << other.coeffs_[ii] << std::endl;
      const_ops_.emplace_back(other.const_ops_[ii]);
      coeffs_.emplace_back(other.coeffs_[ii]);
    }
    return *this;
  } // ... operator+=(...)

  ThisType& operator-=(const BaseType& other)
  {
    LOG_(debug) << "operator-=(other_op=" << &other << ")" << std::endl;
    this->add(other, -1.);
    return *this;
  }

  ThisType& operator-=(const ThisType& other)
  {
    this->logger.state_or(other.logger.state);
    LOG_(debug) << "operator-=(other_const_lincomb_op=" << &other << ")" << std::endl;
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      LOG_(debug) << "  adding op=" << &(other.const_ops_[ii]) << ", coeff=" << -1 * other.coeffs_[ii] << std::endl;
      const_ops_.emplace_back(other.const_ops_[ii]);
      coeffs_.emplace_back(-1 * other.coeffs_[ii]);
    }
    return *this;
  } // ... operator-=(...)

  /// \}
  /// \name These numeric operators override the const ones from OperatorInterface to avoid segfaults due to
  ///       temporaries (which we achieve since *this has the correct type here to select the correct add)
  /// \{

  ConstLincombOperatorType operator*(const FieldType& alpha) const final
  {
    return BaseType::make_operator_mul(*this, alpha);
  }

  ConstLincombOperatorType operator/(const FieldType& alpha) const final
  {
    return BaseType::make_operator_div(*this, alpha);
  }

  ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const final
  {
    return BaseType::make_operator_addsub(*this, other, /*add=*/true);
  }

  ConstLincombOperatorType operator+(const BaseType& other) const final
  {
    return BaseType::make_operator_addsub(*this, other, /*add=*/true);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  ConstLincombOperatorType operator+(const VectorType& vector) const final
  {
    return BaseType::make_operator_addsub(*this, vector, /*add=*/true);
  }

  ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const final
  {
    return BaseType::make_operator_addsub(*this, other, /*add=*/false);
  }

  ConstLincombOperatorType operator-(const BaseType& other) const final
  {
    return BaseType::make_operator_addsub(*this, other, /*add=*/false);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  // we need to implement this to have *this the correct type in the add() below
  ConstLincombOperatorType operator-(const VectorType& vector) const final
  {
    return BaseType::make_operator_addsub(*this, vector, /*add=*/false);
  }

  /// \}

protected:
  const AssemblyGridViewType& assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::vector<std::shared_ptr<OperatorType>> keep_alive_;
  std::vector<XT::Common::ConstStorageProvider<OperatorType>> const_ops_;
  std::vector<FieldType> coeffs_;
}; // class ConstLincombOperator


template <class MatrixType, // <- needs to be manually specified
          class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC>
auto make_const_lincomb_operator(const AssemblyGridViewType& assembly_grid_view,
                                 const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                                 const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                                 const std::string& logging_prefix = "",
                                 const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return ConstLincombOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}

template <class AssemblyGridViewType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
auto make_const_lincomb_operator(const AssemblyGridViewType& assembly_grid_view,
                                 const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                                 const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                                 const std::string& logging_prefix = "",
                                 const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_const_lincomb_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}


template <class MatrixType, // <- needs to be manually specified
          class GV,
          size_t r,
          size_t rC,
          class F>
auto make_const_lincomb_operator(const SpaceInterface<GV, r, r, F>& space,
                                 const std::string& logging_prefix = "",
                                 const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return ConstLincombOperator<GV, r, rC, r, rC, F, MatrixType, GV, GV>(
      space.grid_view(), space, space, logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F>
auto make_const_lincomb_operator(const SpaceInterface<GV, r, r, F>& space,
                                 const std::string& logging_prefix = "",
                                 const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_const_lincomb_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(space, logging_prefix, logging_state);
}


template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AGV,
          class RGV = AGV>
class LincombOperator : public ConstLincombOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
{
  using ThisType = LincombOperator;
  using BaseType = ConstLincombOperator<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

public:
  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::FieldType;
  using typename BaseType::LincombOperatorType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::OperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  LincombOperator(const AssemblyGridViewType& assembly_grid_vw,
                  const SourceSpaceType& src_space,
                  const RangeSpaceType& rng_space,
                  const std::string& logging_prefix = "",
                  const std::array<bool, 3>& logging_enabled = XT::Common::default_logger_state())
    : BaseType(assembly_grid_vw,
               src_space,
               rng_space,
               logging_prefix.empty() ? "LincombOperator" : logging_prefix,
               logging_enabled)
  {
    LOG_(debug) << "LincombOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &src_space
                << ", range_space=" << &rng_space << ")" << std::endl;
  }

  LincombOperator(ThisType& other) = default;

  LincombOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;
  using BaseType::jacobian;
  using BaseType::operator+;
  using BaseType::operator+=;
  using BaseType::operator-;
  using BaseType::operator-=;
  using BaseType::operator*;
  using BaseType::operator*=;
  using BaseType::operator/;
  using BaseType::operator/=;
  using BaseType::add;
  using BaseType::op;

  /// \name Required by BilinearFormInterface.
  /// \{

  void assemble(const bool use_tbb = false) final
  {
    for (auto& oo : ops_)
      oo.access().assemble(use_tbb);
  }

  /// \}
  /// \name These methods extend the ones from ConstLincombOperator
  /// \{

  OperatorType& op(const size_t ii)
  {
    DUNE_THROW_IF(ii >= this->num_ops(),
                  Exceptions::operator_error,
                  "ii = " << ii << "\n   this->num_ops() = " << this->num_ops());
    return ops_[ii].access();
  }

  void add(OperatorType& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op.logger.state);
    LOG_(debug) << "add(op_ref=" << &op << ", coeff=" << coeff << ")" << std::endl;
    ops_.emplace_back(op);
    BaseType::add(ops_.back().access(), coeff); // extends parameter type
  }

  void add(OperatorType*&& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op->logger.state);
    LOG_(debug) << "add(op_ptr=" << op << ", coeff=" << coeff << ")" << std::endl;
    BaseType::add(std::move(op), coeff); // extends parameter type
    ops_.emplace_back(*this->keep_alive_.back());
  }

  void add(ThisType& op, const FieldType& coeff = 1.)
  {
    this->logger.state_or(op.logger.state);
    LOG_(debug) << "add(lincomb_op_ref=" << &op << ", coeff=" << coeff << ")" << std::endl;
    BaseType::add(op, coeff); // // extends parameter type
    for (size_t ii = 0; ii < op.num_ops(); ++ii) {
      LOG_(debug) << "  adding op=" << &(op.ops_[ii]) << std::endl;
      ops_.emplace_back(op.ops_[ii]);
    }
  } // ... add(...)

  // we need this, otherwise add(OperatorType*&&) would be used
  void add(ThisType*&& op, const FieldType& coeff = 1.)
  {
    this->add(*op, coeff);
  }

  /// \}
  /// \name These numeric operators override the ones from ConstLincombOperator to keep track of the mutable ops
  /// \{

  ThisType& operator+=(ThisType& other)
  {
    this->logger.state_or(op.logger.state);
    LOG_(debug) << "operator+=(other_lincomb_op=" << &other << ")" << std::endl;
    this->extend_parameter_type(other.parameter_type());
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      LOG_(debug) << "  adding const_op=" << &(other.const_ops_[ii]) << ", op=" << &(other.ops_[ii])
                  << ", coeff=" << other.coeffs_[ii] << std::endl;
      this->const_ops_.emplace_back(other.const_ops_[ii]);
      ops_.emplace_back(other.ops_[ii]);
      this->coeffs_.emplace_back(other.coeffs_[ii]);
    }
    return *this;
  } // ... operator+=(...)

  ThisType& operator-=(ThisType& other)
  {
    this->logger.state_or(op.logger.state);
    this->extend_parameter_type(other.parameter_type());
    LOG_(debug) << "operator-=(other_lincomb_op=" << &other << ")" << std::endl;
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      LOG_(debug) << "  adding const_op=" << &(other.const_ops_[ii]) << ", op=" << &(other.ops_[ii])
                  << ", coeff=" << -1 * other.coeffs_[ii] << std::endl;
      this->const_ops_.emplace_back(other.const_ops_[ii]);
      ops_.emplace_back(other.ops_[ii]);
      this->coeffs_.emplace_back(-1 * other.coeffs_[ii]);
    }
    return *this;
  } // ... operator-=(...)

  /// \}
  /// \name These numeric operators override the mutable ones from OperatorInterface to avoid segfaults due to
  ///       temporaries (which we achieve since *this has the correct type here to select the correct add)
  /// \{

  LincombOperatorType operator*(const FieldType& alpha)final
  {
    return OperatorType::make_operator_mul(*this, alpha);
  }

  LincombOperatorType operator/(const FieldType& alpha) final
  {
    return OperatorType::make_operator_div(*this, alpha);
  }

  LincombOperatorType operator+(LincombOperatorType& other) final
  {
    return OperatorType::make_operator_addsub(*this, other, /*add=*/true);
  }

  LincombOperatorType operator+(OperatorType& other) final
  {
    return OperatorType::make_operator_addsub(*this, other, /*add=*/true);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  LincombOperatorType operator+(const VectorType& vector) final
  {
    return OperatorType::make_operator_addsub(*this, vector, /*add=*/true);
  }

  LincombOperatorType operator-(LincombOperatorType& other) final
  {
    return OperatorType::make_operator_addsub(*this, other, /*add=*/false);
  }

  LincombOperatorType operator-(OperatorType& other) final
  {
    return OperatorType::make_operator_addsub(*this, other, /*add=*/false);
  }

  /// \note vector is interpreted as a ConstantOperator
  /// \sa ConstantOperator
  LincombOperatorType operator-(const VectorType& vector) final
  {
    return OperatorType::make_operator_addsub(*this, vector, /*add=*/false);
  }

  /// \}

private:
  std::vector<XT::Common::StorageProvider<OperatorType>> ops_;
}; // class LincombOperator


template <class MatrixType, // <- needs to be manually specified
          class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC>
auto make_lincomb_operator(const AssemblyGridViewType& assembly_grid_view,
                           const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                           const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return LincombOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, MatrixType, SGV, RGV>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}

template <class AssemblyGridViewType, class SGV, size_t s_r, size_t s_rC, class F, class RGV, size_t r_r, size_t r_rC>
auto make_lincomb_operator(const AssemblyGridViewType& assembly_grid_view,
                           const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                           const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_lincomb_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(
      assembly_grid_view, source_space, range_space, logging_prefix, logging_state);
}


template <class MatrixType, // <- needs to be manually specified
          class GV,
          size_t r,
          size_t rC,
          class F>
auto make_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return LincombOperator<GV, r, rC, r, rC, F, MatrixType, GV, GV>(
      space.grid_view(), space, space, logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F>
auto make_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space,
                           const std::string& logging_prefix = "",
                           const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  return make_lincomb_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(space, logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LINCOMB_HH
