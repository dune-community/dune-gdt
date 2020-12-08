// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)

#ifndef DUNE_GDT_OPERATORS_CONSTANT_HH
#define DUNE_GDT_OPERATORS_CONSTANT_HH

#include <dune/xt/common/memory.hh>
#include <dune/xt/la/container.hh>
#include <dune/xt/la/type_traits.hh>

#include <dune/gdt/exceptions.hh>
#include <dune/gdt/print.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/// \todo add zero_threshhold?
template <class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class V = XT::LA::IstlDenseVector<F>,
          class RGV = SGV>
class ConstantForwardOperator : public ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, V, RGV>
{
public:
  using ThisType = ConstantForwardOperator;
  using BaseType = ForwardOperatorInterface<SGV, s_r, s_rC, r_r, r_rC, F, V, RGV>;

  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::VectorType;

  ConstantForwardOperator(const RangeSpaceType& rng_space,
                          const VectorType& val,
                          const std::string& logging_prefix = "",
                          const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({}, logging_prefix.empty() ? "ConstantForwardOperator" : logging_prefix, logging_state)
    , range_space_(rng_space)
    , value_(val)
    , is_zero_(XT::Common::is_zero(value_.access().sup_norm()))
  {
    LOG_(debug) << "ConstantForwardOperator(range_space=" << &rng_space
                << ", value.sup_norm()=" << value_.access().sup_norm() << ")" << std::endl;
    this->assert_matching_range(value_.access());
  }

  ConstantForwardOperator(const RangeSpaceType& rng_space,
                          VectorType*&& val,
                          const std::string& logging_prefix = "",
                          const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({}, logging_prefix.empty() ? "ConstantForwardOperator" : logging_prefix, logging_state)
    , range_space_(rng_space)
    , value_(std::move(val))
    , is_zero_(XT::Common::is_zero(value_.access().sup_norm()))
  {
    LOG_(debug) << "ConstantForwardOperator(range_space=" << &rng_space
                << ", value.sup_norm()=" << value_.access().sup_norm() << ")" << std::endl;
    this->assert_matching_range(value_.access());
  }

  ConstantForwardOperator(const ThisType& other) = default;

  ConstantForwardOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  bool linear() const override final
  {
    return is_zero_;
  }

  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    LOG_(info) << "setting range_vector to " << (is_zero_ ? "0" : "constant value") << " ..." << std::endl;
    range_vector = value_.access();
  } // ... apply(...)

  /// \}

private:
  const RangeSpaceType& range_space_;
  const XT::Common::ConstStorageProvider<VectorType> value_;
  const bool is_zero_;
}; // class ConstantForwardOperator


template <class GV, size_t r, size_t rC, class F, class V>
auto make_constant_forward_operator(const SpaceInterface<GV, r, rC, F>& space,
                                    const XT::LA::VectorInterface<V>& value,
                                    const std::string& logging_prefix = "",
                                    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  using V_ = typename XT::LA::VectorInterface<V>::derived_type;
  return ConstantForwardOperator<GV, r, rC, r, rC, F, V_, GV>(space, value.as_imp(), logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F, class VectorType>
auto make_constant_forward_operator(const SpaceInterface<GV, r, rC, F>& space,
                                    VectorType*&& value_ptr,
                                    const std::string& logging_prefix = "",
                                    const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  return ConstantForwardOperator<GV, r, rC, r, rC, F, VectorType, GV>(
      space, std::move(value_ptr), logging_prefix, logging_state);
}


/// \todo add zero_threshhold?
/// \todo only set those dofs to zero which are associated with assembly_grid_view?
template <class AGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class F = double,
          class M = XT::LA::IstlRowMajorSparseMatrix<F>,
          class SGV = AGV,
          class RGV = AGV>
class ConstantOperator : public OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>
{
  using ThisType = ConstantOperator;
  using BaseType = OperatorInterface<AGV, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>;

public:
  using typename BaseType::AssemblyGridViewType;
  using typename BaseType::RangeSpaceType;
  using typename BaseType::SourceFunctionType;
  using typename BaseType::SourceSpaceType;
  using typename BaseType::VectorType;

  ConstantOperator(const AssemblyGridViewType& assembly_grid_vw,
                   const SourceSpaceType& src_space,
                   const RangeSpaceType& rng_space,
                   const VectorType& val,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({}, logging_prefix.empty() ? "ConstantOperator" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_vw)
    , source_space_(src_space)
    , range_space_(rng_space)
    , value_(val)
    , is_zero_(XT::Common::is_zero(value_.access().sup_norm()))
  {
    LOG_(debug) << "ConstantOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &src_space
                << ", range_space=" << &rng_space << ", value.sup_norm()=" << value_.access().sup_norm() << ")"
                << std::endl;
    this->assert_matching_range(value_.access());
  }

  ConstantOperator(const AssemblyGridViewType& assembly_grid_vw,
                   const SourceSpaceType& src_space,
                   const RangeSpaceType& rng_space,
                   VectorType*&& val,
                   const std::string& logging_prefix = "",
                   const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
    : BaseType({}, logging_prefix.empty() ? "ConstantOperator" : logging_prefix, logging_state)
    , assembly_grid_view_(assembly_grid_vw)
    , source_space_(src_space)
    , range_space_(rng_space)
    , value_(std::move(val))
    , is_zero_(XT::Common::is_zero(value_.access().sup_norm()))
  {
    LOG_(debug) << "ConstantOperator(assembly_grid_view=" << &assembly_grid_vw << ", source_space=" << &src_space
                << ", range_space=" << &rng_space << ", value.sup_norm()=" << value_.access().sup_norm() << ")"
                << std::endl;
    this->assert_matching_range(value_.access());
  }

  ConstantOperator(const ThisType& other) = default;

  ConstantOperator(ThisType&& source) = default;

  // pull in methods from various base classes
  using BaseType::apply;

  /// \name Required by ForwardOperatorInterface.
  /// \{

  const RangeSpaceType& range_space() const override final
  {
    return range_space_;
  }

  bool linear() const override final
  {
    return is_zero_;
  }

  // avoid non-optimal default implementation in OperatorInterface
  void apply(SourceFunctionType source_function,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "apply(source_function=" << &source_function
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_range(range_vector);
    LOG_(info) << "setting range_vector to " << (is_zero_ ? "0" : "constant value") << " ..." << std::endl;
    range_vector = value_.access();
  } // ... apply(...)

  /// \}
  /// \name Required by OperatorInterface.
  /// \{

  const SourceSpaceType& source_space() const override final
  {
    return source_space_;
  }

  const AssemblyGridViewType& assembly_grid_view() const override final
  {
    return assembly_grid_view_;
  }

  void apply(const VectorType& source_vector,
             VectorType& range_vector,
             const XT::Common::Parameter& param = {}) const override final
  {
    LOG_(debug) << "apply(source_vector.sup_norm()=" << source_vector.sup_norm()
                << ", range_vector.sup_norm()=" << range_vector.sup_norm() << ", param=" << print(param) << ")"
                << std::endl;
    this->assert_matching_source(source_vector);
    this->assert_matching_range(range_vector);
    LOG_(info) << "setting range_vector to " << (is_zero_ ? "0" : "constant value") << " ..." << std::endl;
    range_vector = value_.access();
  } // ... apply(...)

protected:
  std::vector<XT::Common::Configuration> all_jacobian_options() const override final
  {
    return {{{"type", "zero"}}};
  }

  /// \}

private:
  const AssemblyGridViewType& assembly_grid_view_;
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  const XT::Common::ConstStorageProvider<VectorType> value_;
  const bool is_zero_;
}; // class ConstantOperator


template <class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class V>
auto make_constant_operator(const AssemblyGridViewType& assembly_grid_view,
                            const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                            const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                            const XT::LA::VectorInterface<V>& value,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  using M = XT::LA::matrix_t<typename XT::LA::VectorInterface<V>::derived_type>;
  return ConstantOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>(
      assembly_grid_view, source_space, range_space, value.as_imp(), logging_prefix, logging_state);
}

template <class AssemblyGridViewType,
          class SGV,
          size_t s_r,
          size_t s_rC,
          class F,
          class RGV,
          size_t r_r,
          size_t r_rC,
          class VectorType>
auto make_constant_operator(const AssemblyGridViewType& assembly_grid_view,
                            const SpaceInterface<SGV, s_r, s_rC, F>& source_space,
                            const SpaceInterface<RGV, r_r, r_rC, F>& range_space,
                            const VectorType*&& value_ptr,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::Grid::is_view<AssemblyGridViewType>::value, "");
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  using M = XT::LA::matrix_t<VectorType>;
  return ConstantOperator<AssemblyGridViewType, s_r, s_rC, r_r, r_rC, F, M, SGV, RGV>(
      assembly_grid_view, source_space, range_space, std::move(value_ptr), logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F, class V>
auto make_constant_operator(const SpaceInterface<GV, r, rC, F>& space,
                            const XT::LA::VectorInterface<V>& value,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  using M = XT::LA::matrix_t<typename XT::LA::VectorInterface<V>::derived_type>;
  return ConstantOperator<GV, r, rC, r, rC, F, M, GV, GV>(
      space.grid_view(), space, space, value.as_imp(), logging_prefix, logging_state);
}

template <class GV, size_t r, size_t rC, class F, class VectorType>
auto make_constant_operator(const SpaceInterface<GV, r, rC, F>& space,
                            VectorType*&& value_ptr,
                            const std::string& logging_prefix = "",
                            const std::array<bool, 3>& logging_state = XT::Common::default_logger_state())
{
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  using M = XT::LA::matrix_t<VectorType>;
  return ConstantOperator<GV, r, rC, r, rC, F, M, GV, GV>(
      space.grid_view(), space, space, std::move(value_ptr), logging_prefix, logging_state);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_CONSTANT_HH
