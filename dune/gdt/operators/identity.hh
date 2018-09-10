// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reseVed.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_OPERATORS_IDENTITY_HH
#define DUNE_GDT_OPERATORS_IDENTITY_HH

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/bilinear-forms/generic.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * See also OperatorInterface for a description of the template arguemnts.
 *
 * \sa OperatorInterface
 */
template <class M, class GV, size_t r = 1, size_t rC = 1>
class IdentityOperator : public OperatorInterface<M, GV, r, rC>
{
  using ThisType = IdentityOperator<M, GV, r, rC>;
  using BaseType = OperatorInterface<M, GV, r, rC>;

public:
  using typename BaseType::SourceSpaceType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::VectorType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::F;

  IdentityOperator(const SourceSpaceType& spc)
    : space_(spc)
  {
  }

  bool linear() const override final
  {
    return true;
  }

  const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  using BaseType::apply;

  void
  apply(const VectorType& source, VectorType& range, const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!space_.contains(source), Exceptions::operator_error, "");
    DUNE_THROW_IF(!space_.contains(range), Exceptions::operator_error, "");
    range = source;
  }

  std::vector<std::string> invert_options() const override final
  {
    return {"identity"};
  }

  XT::Common::Configuration invert_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->invert_options().at(0), Exceptions::operator_error, "type = " << type);
    return {{"type", type}};
  }

  using BaseType::apply_inverse;

  void apply_inverse(const VectorType& range,
                     VectorType& source,
                     const XT::Common::Configuration& opts,
                     const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!space_.contains(range), Exceptions::operator_error, "");
    DUNE_THROW_IF(!space_.contains(source), Exceptions::operator_error, "");
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, "");
    DUNE_THROW_IF(
        opts.get<std::string>("type") != this->invert_options().at(0), Exceptions::operator_error, "opts = " << opts);
    source = range;
  }

  std::vector<std::string> jacobian_options() const override final
  {
    return {"identity"};
  }

  XT::Common::Configuration jacobian_options(const std::string& type) const override final
  {
    DUNE_THROW_IF(type != this->jacobian_options().at(0), Exceptions::operator_error, "type = " << type);
    return {{"type", type}};
  }

  using BaseType::jacobian;

  void jacobian(const VectorType& /*source*/,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(jacobian_op.source_space().mapper().size() != space_.mapper().size(), Exceptions::operator_error, "");
    DUNE_THROW_IF(jacobian_op.range_space().mapper().size() != space_.mapper().size(), Exceptions::operator_error, "");
    DUNE_THROW_IF(!opts.has_key("type"), Exceptions::operator_error, "opts = \n");
    DUNE_THROW_IF(
        opts.get<std::string>("type") != this->jacobian_options().at(0), Exceptions::operator_error, "opts = " << opts);
    jacobian_op.append(GenericLocalElementBilinearForm<XT::Grid::extract_entity_t<GV>, r, rC, F>(
        [](const auto& test_basis, const auto& /*ansatz_basis*/, auto& result, const auto& param) {
          for (size_t ii = 0; ii < test_basis.size(param); ++ii)
            result[ii][ii] = 1.;
        }));
  } // ... jacobian(...)

private:
  const SourceSpaceType& space_;
}; // class IdentityOperator


template <class Matrix, class GV, size_t r, size_t rC, class F>
std::enable_if_t<XT::LA::is_matrix<Matrix>::value, IdentityOperator<Matrix, GV, r, rC>>
make_identity_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return IdentityOperator<Matrix, GV, r, rC>(space);
}


template <class GV, size_t r, size_t rC, class F>
IdentityOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>
make_identity_operator(const SpaceInterface<GV, r, rC, F>& space)
{
  return IdentityOperator<typename XT::LA::Container<F>::MatrixType, GV, r, rC>(space);
}


template <class M, class GV, size_t r, size_t rC>
IdentityOperator<M, GV, r, rC> make_identity_operator(const OperatorInterface<M, GV, r, rC, r, rC, GV>& op)
{
  return IdentityOperator<M, GV, r, rC>(op.source_space());
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_IDENTITY_HH
