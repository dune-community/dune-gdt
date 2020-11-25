// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_IDENTITY_HH
#define DUNE_GDT_OPERATORS_IDENTITY_HH

#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/local/operators/generic.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


template <class GV, size_t r = 1, size_t rC = 1, class F = double, class M = XT::LA::IstlRowMajorSparseMatrix<F>>
class IdentityOperator : public OperatorInterface<GV, r, rC, r, rC, F, M, GV, GV>
{
public:
  using ThisType = IdentityOperator;
  using BaseType = OperatorInterface<GV, r, rC, r, rC, F, M, GV, GV>;

  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::MatrixOperatorType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  IdentityOperator(const SourceSpaceType& space,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = {{false, false, true}})
    : BaseType({}, logging_prefix.empty() ? "IdentityOperator" : logging_prefix, logging_state)
    , space_(space)
  {
    LOG_(debug) << "IdentityOperator(space=" << &space << ")" << std::endl;
  }

  // pull in methods from various base classes
  using BaseType::apply;
  using BaseType::apply_inverse;
  using BaseType::jacobian;

  /// \name Required by OperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return space_;
  }

  bool linear() const override final
  {
    return true;
  }

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return space_;
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return space_.grid_view();
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "apply(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << param << ")" << std::endl;
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    LOG_(info) << "setting range_vector = source_vector ..." << std::endl;
    range_vector = source_vector;
  } // ... apply(...)

protected:
  std::vector<XT::Common::Configuration> all_jacobian_options() const override final
  {
    return {{{"type", "identity"}}};
  }

  virtual std::vector<XT::Common::Configuration> all_invert_options() const
  {
    return {{{"type", "identity"}}};
  }

public:
  void jacobian(const VectorType& source_vector,
                MatrixOperatorType& jacobian_op,
                const XT::Common::Configuration& opts,
                const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "jacobian(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", jacobian_op.matrix().sup_norm()=" << jacobian_op.matrix().sup_norm()
                << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << param << ")" << std::endl;
    this->assert_jacobian_opts(opts); // ensures that type identity is requested
    LOG_(info) << "adding unit diagonal to jacobian_op ..." << std::endl;
    const F unit = jacobian_op.scaling;
    const size_t size = jacobian_op.matrix().rows();
    DUNE_THROW_IF(jacobian_op.matrix().cols() != size, Exceptions::operator_error, "jacobian_op is not square!");
    for (size_t ii = 0; ii < size; ++ii)
      jacobian_op.matrix().add_to_entry(ii, ii, unit);
  } // ... jacobian(...)

  void apply_inverse(const VectorType& range_vector,
                     VectorType& source_vector,
                     const XT::Common::Configuration& opts,
                     const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "apply_inverse(range_vector.sup_norm()=" << source_vector.sup_norm()
                << ", source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", opts=" << print(opts, {{"oneline", "true"}}) << ", param=" << param << ")" << std::endl;
    this->assert_apply_inverse_opts(opts); // ensures that type identity is requested
    LOG_(info) << "setting source_vector = range_vector ..." << std::endl;
    source_vector = range_vector;
  } // ... apply_inverse(...)

  /// \}

private:
  const SourceSpaceType& space_;
}; // namespace GDT


template <class MatrixType, // <- needs to be manually specified
          class GV,
          size_t r,
          size_t rC,
          class F>
auto make_identity_operator(const SpaceInterface<GV, r, rC, F>& space,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  return IdentityOperator<GV, r, rC, F, MatrixType>(space.grid_view(), space, space, logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F>
auto make_identity_operator(const SpaceInterface<GV, r, rC, F>& space,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  return make_identity_operator<XT::LA::IstlRowMajorSparseMatrix<F>>(space, logging_prefix, logging_state);
}


template <class GV, size_t r, size_t rC, class F, class M>
auto make_identity_operator(const OperatorInterface<GV, r, rC, r, rC, F, M, GV, GV>& some_operator,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = {{false, false, true}})
{
  // check if source and range coincide
  const auto& source = some_operator.source_space();
  const auto& range = some_operator.range_space();
  DUNE_THROW_IF(source.type() != range.type(),
                Exceptions::operator_error,
                "Can not create IdentityOperator like given operator, source and range do not coincide!"
                    << "\n"
                    << "   source.type() = " << source.type() << "\n"
                    << "   range.type() = " << range.type());
  DUNE_THROW_IF(source.mapper().size() != range.mapper().size(),
                Exceptions::operator_error,
                "Can not create IdentityOperator like given operator, source and range do not coincide!"
                    << "\n"
                    << "   source.mapper().size() = " << source.mapper().size() << "\n"
                    << "   range.mapper().size() = " << range.mapper().size());
  // we still cannot rule out different source and range, but have no further means to check
  return IdentityOperator<GV, r, rC, F, M>(source, logging_prefix, logging_state);
} // ... make_identity_operator(...)


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_IDENTITY_HH
