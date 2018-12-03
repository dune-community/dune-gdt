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


/**
 * \todo Add handling of parametric operators: merge this->parameter_type() and op.parameter_type() upon add().
 */
template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class ConstLincombOperator : public OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
{
  using ThisType = ConstLincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using BaseType = OperatorInterface<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

public:
  using typename BaseType::ConstLincombOperatorType;
  using typename BaseType::FieldType;
  using typename BaseType::LincombOperatorType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  using OperatorType = BaseType;

  ConstLincombOperator(const SourceSpaceType& src_space, const RangeSpaceType& rng_space)
    : source_space_(src_space)
    , range_space_(rng_space)
  {}

  ConstLincombOperator(const ThisType& other) = default;

  ConstLincombOperator(ThisType&& source) = default;

  void add(const OperatorType& op, const FieldType& coeff = 1.)
  {
    const_ops_.emplace_back(op);
    coeffs_.emplace_back(coeff);
  }

  void add(OperatorType*&& op, const FieldType& coeff = 1.)
  {
    keep_alive_.emplace_back(std::move(op));
    const_ops_.emplace_back(*keep_alive_.back());
    coeffs_.emplace_back(coeff);
  }

  void add(const ThisType& op, const FieldType& coeff = 1.)
  {
    // Only adding op itself would lead to segfaults if op is a temporary
    for (size_t ii = 0; ii < op.num_ops(); ++ii) {
      const_ops_.emplace_back(op.const_ops_[ii]);
      coeffs_.emplace_back(coeff * op.coeffs_[ii]);
    }
  }

  void add(ThisType*&& op, const FieldType& coeff = 1.)
  {
    this->add(*op, coeff);
  }

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

  bool linear() const override final
  {
    for (const auto& op : const_ops_)
      if (!op.access().linear())
        return false;
    return true;
  }

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  using BaseType::apply;

  void apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& param = {}) const override final
  {
    range.set_all(0);
    auto tmp = range;
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      this->op(ii).apply(source, tmp, param);
      tmp *= this->coeff(ii);
      range += tmp;
    }
  } // ... append(...)

  std::vector<std::string> jacobian_options() const override final
  {
    return {"lincomb"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->jacobian_options().at(0), Exceptions::operator_error, "type = " << type);
    using XT::Common::to_string;
    XT::Common::Configuration ret({{"type", type}});
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      for (auto&& tp : this->op(ii).jacobian_options())
        ret["op_" + to_string(ii) + ".types"] += ", " + tp;
      ret.add(this->op(ii).jacobian_options(this->op(ii).jacobian_options().at(0)), "op_" + to_string(ii) + ".opts");
    }
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      const auto back_ii = ssize_t(this->num_ops()) - 1 - ssize_t(ii);
      for (auto&& tp : this->op(back_ii).jacobian_options())
        ret["back_op_" + to_string(back_ii) + ".types"] += ", " + tp;
      ret.add(this->op(back_ii).jacobian_options(this->op(back_ii).jacobian_options().at(0)),
              "back_op_" + to_string(back_ii) + ".opts");
    }
    return ret;
  } // ... jacobian_options(...)

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
    using XT::Common::to_string;
    const XT::Common::Configuration default_opts = this->jacobian_options(this->jacobian_options().at(0));
    for (size_t ii = 0; ii < this->num_ops(); ++ii) {
      // parse opts
      XT::Common::Configuration op_opts;
      if (opts.has_sub("op_" + to_string(ii) + ".opts")) {
        if (opts.sub("op_" + to_string(ii) + ".opts") != default_opts.sub("op_" + to_string(ii) + ".opts"))
          op_opts = opts.sub("op_" + to_string(ii) + ".opts");
      }
      const auto back_ii = ssize_t(this->num_ops()) - 1 - ssize_t(ii);
      if (opts.has_sub("back_op_" + to_string(back_ii) + ".opts")) {
        if (opts.sub("back_op_" + to_string(back_ii) + ".opts")
            != default_opts.sub("back_op_" + to_string(back_ii) + ".opts")) {
          DUNE_THROW_IF(opts.has_sub("op_" + to_string(ii) + ".opts"),
                        Exceptions::operator_error,
                        "Cannot define opts for operator " << ii << " by specifying op_" << ii << ".opts and back_op_"
                                                           << (back_ii) << ".opts at the same time!\n\n\nopts = \n"
                                                           << opts);
          op_opts = opts.sub("back_op_" + to_string(back_ii) + ".opts");
        }
      }
      if (op_opts.empty())
        op_opts["type"] = this->op(ii).jacobian_options().at(0);
      // save curent scaling
      const auto scaling = jacobian_op.scaling;
      // add op
      jacobian_op.scaling *= this->coeff(ii);
      this->op(ii).jacobian(source, jacobian_op, op_opts, param);
      // restore scaling
      jacobian_op.scaling = scaling;
    }
  } // ... jacobian(...)

  ThisType& operator*=(const FieldType& alpha)
  {
    for (auto& coeff : coeffs_)
      coeff *= alpha;
    return *this;
  }

  ThisType& operator/=(const FieldType& alpha)
  {
    for (auto& coeff : coeffs_)
      coeff /= alpha;
    return *this;
  }

  ThisType& operator+=(const BaseType& other)
  {
    this->add(other);
    return *this;
  }

  ThisType& operator+=(const ThisType& other)
  {
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      const_ops_.emplace_back(other.const_ops_[ii]);
      coeffs_.emplace_back(other.coeffs_[ii]);
    }
    return *this;
  }

  ThisType& operator-=(const BaseType& other)
  {
    this->add(other, -1.);
    return *this;
  }

  ThisType& operator-=(const ThisType& other)
  {
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      const_ops_.emplace_back(other.const_ops_[ii]);
      coeffs_.emplace_back(-1 * other.coeffs_[ii]);
    }
    return *this;
  }

  // We need to override some operator+-*/ from the interface to avoid segfaults due to temporaries

  ConstLincombOperatorType operator*(const FieldType& alpha) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret *= alpha;
    return ret;
  }

  ConstLincombOperatorType operator/(const FieldType& alpha) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret /= alpha;
    return ret;
  }

  using BaseType::operator+;

  ConstLincombOperatorType operator+(const ConstLincombOperatorType& other) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  ConstLincombOperatorType operator+(const BaseType& other) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  ConstLincombOperatorType operator+(const VectorType& vector) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            1.);
    return ret;
  }

  using BaseType::operator-;

  ConstLincombOperatorType operator-(const ConstLincombOperatorType& other) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  ConstLincombOperatorType operator-(const BaseType& other) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  ConstLincombOperatorType operator-(const VectorType& vector) const override final
  {
    ConstLincombOperatorType ret(*this);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            -1.);
    return ret;
  }

protected:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  std::vector<std::shared_ptr<OperatorType>> keep_alive_;
  std::vector<XT::Common::ConstStorageProvider<OperatorType>> const_ops_;
  std::vector<FieldType> coeffs_;
}; // class ConstLincombOperator


template <class Matrix, class GV, size_t r, size_t rC, class F>
std::enable_if_t<XT::LA::is_matrix<Matrix>::value, ConstLincombOperator<Matrix, GV, r, rC>>
make_const_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return ConstLincombOperator<Matrix, GV, r, rC>(space, space);
}


template <class GV, size_t r, size_t rC, class F>
ConstLincombOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>
make_const_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return ConstLincombOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>(space, space);
}


template <class M, class SGV, size_t s_r = 1, size_t s_rC = 1, size_t r_r = s_r, size_t r_rC = s_rC, class RGV = SGV>
class LincombOperator : public ConstLincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>
{
  using ThisType = LincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;
  using BaseType = ConstLincombOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>;

public:
  using typename BaseType::FieldType;
  using typename BaseType::LincombOperatorType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::OperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  LincombOperator(const SourceSpaceType& src_space, const RangeSpaceType& rng_space)
    : BaseType(src_space, rng_space)
  {}

  LincombOperator(ThisType& other)
    : BaseType(other)
  {
    for (auto& op : other.ops_)
      this->ops_.emplace_back(op);
  }

  LincombOperator(ThisType&& source)
    : BaseType(source)
    , ops_(std::move(source.ops_))
  {}

  using BaseType::add;

  void add(OperatorType& op, const FieldType& coeff = 1.)
  {
    ops_.emplace_back(op);
    BaseType::add(ops_.back().access(), coeff);
  }

  void add(OperatorType*&& op, const FieldType& coeff = 1.)
  {
    BaseType::add(std::move(op), coeff);
    ops_.emplace_back(*this->keep_alive_.back());
  }

  void add(ThisType& op, const FieldType& coeff = 1.)
  {
    BaseType::add(op, coeff);
    for (size_t ii = 0; ii < op.num_ops(); ++ii)
      ops_.emplace_back(op.ops_[ii]);
  }

  void add(ThisType*&& op, const FieldType& coeff = 1.)
  {
    this->add(*op, coeff);
  }

  using BaseType::op;

  OperatorType& op(const size_t ii)
  {
    DUNE_THROW_IF(ii >= this->num_ops(),
                  Exceptions::operator_error,
                  "ii = " << ii << "\n   this->num_ops() = " << this->num_ops());
    return ops_[ii].access();
  }

  OperatorType& assemble(const bool use_tbb = false) override final
  {
    for (auto& op : ops_)
      op.access().assemble(use_tbb);
    return *this;
  }

  // we need to override some operators, see above

  using BaseType::operator+=;

  ThisType& operator+=(ThisType& other)
  {
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      this->const_ops_.emplace_back(other.const_ops_[ii]);
      ops_.emplace_back(other.ops_[ii]);
      this->coeffs_.emplace_back(other.coeffs_[ii]);
    }
    return *this;
  }

  using BaseType::operator-=;

  ThisType& operator-=(ThisType& other)
  {
    for (size_t ii = 0; ii < other.num_ops(); ++ii) {
      this->const_ops_.emplace_back(other.const_ops_[ii]);
      ops_.emplace_back(other.ops_[ii]);
      this->coeffs_.emplace_back(-1 * other.coeffs_[ii]);
    }
    return *this;
  }

  using BaseType::operator*;

  LincombOperatorType operator*(const FieldType& alpha)override final
  {
    LincombOperatorType ret(*this);
    ret *= alpha;
    return ret;
  }

  using BaseType::operator/;

  LincombOperatorType operator/(const FieldType& alpha) override final
  {
    LincombOperatorType ret(*this);
    ret /= alpha;
    return ret;
  }

  using BaseType::operator+;

  LincombOperatorType operator+(LincombOperatorType& other) override final
  {
    LincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  LincombOperatorType operator+(OperatorType& other) override final
  {
    LincombOperatorType ret(*this);
    ret += other;
    return ret;
  }

  LincombOperatorType operator+(const VectorType& vector) override final
  {
    LincombOperatorType ret(*this);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            -1.);
    return ret;
  }

  using BaseType::operator-;

  LincombOperatorType operator-(LincombOperatorType& other) override final
  {
    LincombOperatorType ret(*this);
    ret -= other;
    return ret;
  }

  LincombOperatorType operator-(OperatorType& other) override final
  {
    LincombOperatorType ret(*this);
    ret -= other;
    return ret;
  }

  LincombOperatorType operator-(const VectorType& vector) override final
  {
    LincombOperatorType ret(*this);
    ret.add(new ConstantOperator<M, SGV, s_r, s_rC, r_r, r_rC, RGV>(this->source_space(), this->range_space(), vector),
            -1.);
    return ret;
  }

private:
  std::vector<XT::Common::StorageProvider<OperatorType>> ops_;
}; // class LincombOperator


template <class Matrix, class GV, size_t r, size_t rC, class F>
std::enable_if_t<XT::LA::is_matrix<Matrix>::value, LincombOperator<Matrix, GV, r, rC>>
make_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return LincombOperator<Matrix, GV, r, rC>(space, space);
}


template <class GV, size_t r, size_t rC, class F>
LincombOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>
make_lincomb_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return LincombOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>(space, space);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_LINCOMB_HH
