// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_FD_JACOBIAN_ASSEMBLERS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_FD_JACOBIAN_ASSEMBLERS_HH

#include <cmath>
#include <memory>

#include <dune/common/dynvector.hh>

#include <dune/xt/common/float_cmp.hh>
#include <dune/xt/common/parameter.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/local/operators/interfaces.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Computes a finite-difference approximation of the jacobian of an operator induced by a local element operator.
 *
 * \note Presumes that the nonlinearity in the first argument of the operator does not suffer from restriction to the
 *       neighborhood.
 *
 * \note This implementation is not optimal, since it requires a full source and range vector. This can only be fixed
 *       after refactoring local discrete functions and local dof vectors.
 *
 * See also LocalElementOperatorInterface for a description of the template arguments.
 *
 * \sa LocalElementOperatorInterface
 */
template <class M,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class F = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RGV = SGV>
class LocalElementOperatorFiniteDifferenceJacobianAssembler : public XT::Grid::ElementFunctor<SGV>
{
  static_assert(XT::LA::is_matrix<M>::value, "");
  static_assert(XT::Grid::is_view<SGV>::value, "");
  static_assert(XT::Grid::is_view<RGV>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<SGV>, XT::Grid::extract_entity_t<RGV>>::value, "");

  using BaseType = XT::Grid::ElementFunctor<SGV>;
  using ThisType = LocalElementOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC, RGV>;

public:
  using typename BaseType::ElementType;

  using MatrixType = M;
  using VectorType = XT::LA::vector_t<M>;
  using V = VectorType;

  using LocalElementOperatorType = LocalElementOperatorInterface<V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV>;
  using SourceSpaceType = typename LocalElementOperatorType::DiscreteSourceType::SpaceType;
  using RangeSpaceType = typename LocalElementOperatorType::LocalRangeType::SpaceType;

  LocalElementOperatorFiniteDifferenceJacobianAssembler(const SourceSpaceType& source_space,
                                                        const RangeSpaceType& range_space,
                                                        MatrixType& matrix,
                                                        const VectorType& source_vector,
                                                        const LocalElementOperatorType& local_operator,
                                                        const XT::Common::Parameter& param = {})
    : BaseType()
    , source_space_(source_space)
    , range_space_(range_space)
    , matrix_(matrix)
    , source_vector_(source_vector)
    , param_(param)
    , scaling_(param_.has_key("matrixoperator.scaling") ? param_.get("matrixoperator.scaling").at(0) : 1.)
    , eps_(param_.has_key("finite-difference-jacobians.eps") ? param_.get("finite-difference-jacobians.eps").at(0)
                                                             : 1e-7)
    , source_(source_space_) // This is a full vector, and intended!
    , range_(range_space_) // This is a full vector, and intended!
    , local_source_(source_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , local_op_(local_operator.with_source(source_))
  {
    source_.dofs().vector() = source_vector_;
  }

  LocalElementOperatorFiniteDifferenceJacobianAssembler(const ThisType& other)
    : BaseType(other)
    , source_space_(other.source_space_)
    , range_space_(other.range_space_)
    , matrix_(other.matrix_)
    , source_vector_(other.source_vector_)
    , param_(other.param_)
    , scaling_(other.scaling_)
    , eps_(other.eps_)
    , source_(source_space_) // This is a full vector, and intended!
    , range_(range_space_) // This is a full vector, and intended!
    , local_source_(source_.local_discrete_function())
    , local_range_(range_.local_discrete_function())
    , local_op_(other.local_op_->with_source(source_))
  {
    source_.dofs().vector() = source_vector_;
  }

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const ElementType& element) override final
  {
    // some preparations
    local_source_->bind(element);
    local_range_->bind(element);
    source_space_.mapper().global_indices(element, global_source_indices_);
    range_space_.mapper().global_indices(element, global_range_indices_);
    const size_t local_source_size = source_space_.mapper().local_size(element);
    const size_t local_range_size = range_space_.mapper().local_size(element);
    if (range_DoFs_.size() < local_range_size)
      range_DoFs_.resize(local_range_size, 0);
    local_range_->dofs().set_all(0);
    // apply op as is, keep the result, clear local range
    local_op_->apply(*local_range_, param_);
    for (size_t ii = 0; ii < local_range_size; ++ii)
      range_DoFs_[ii] = local_range_->dofs()[ii];
    local_range_->dofs().set_all(0);
    // loop over all source DoFs
    for (size_t jj = 0; jj < local_source_size; ++jj) {
      // perturb source DoF
      const auto jjth_source_DoF = local_source_->dofs()[jj];
      const auto eps = eps_ * (1. + std::abs(jjth_source_DoF));
      local_source_->dofs()[jj] += eps;
      // apply op with perturbed source DoF
      local_op_->apply(*local_range_, param_);
      // observe perturbation in range DoFs
      for (size_t ii = 0; ii < local_range_size; ++ii) {
        auto derivative = (local_range_->dofs()[ii] - range_DoFs_[ii]) / eps;
        if (XT::Common::FloatCmp::eq(derivative, eps))
          derivative = 0;
        matrix_.add_to_entry(global_range_indices_[ii], global_source_indices_[jj], scaling_ * derivative);
      }
      // restore source
      local_source_->dofs()[jj] = jjth_source_DoF;
    }
  } // ... apply_local(...)

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  MatrixType& matrix_;
  const VectorType& source_vector_;
  const XT::Common::Parameter param_;
  const double scaling_;
  const real_t<F> eps_;
  DiscreteFunction<V, SGV, s_r, s_rC, F> source_;
  DiscreteFunction<V, RGV, r_r, r_rC, F> range_;
  std::unique_ptr<LocalDiscreteFunction<V, SGV, s_r, s_rC, F>> local_source_;
  std::unique_ptr<LocalDiscreteFunction<V, RGV, r_r, r_rC, F>> local_range_;
  DynamicVector<size_t> global_source_indices_;
  DynamicVector<size_t> global_range_indices_;
  DynamicVector<F> range_DoFs_;
  const std::unique_ptr<LocalElementOperatorType> local_op_;
}; // class LocalElementOperatorFiniteDifferenceJacobianAssembler


/**
 * \brief Computes a finite-difference approximation of the jacobian of an operator induced by a local intersection
 *        operator.
 *
 * \note Presumes that the nonlinearity in the first argument of the operator does not suffer from restriction to the
 *       neighborhood.
 *
 * \note This implementation is not optimal, since it requires a full source and range vector. This can only be fixed
 *       after refactoring local discrete functions and local dof vectors.
 *
 * See also LocalIntersectionOperatorInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionOperatorInterface
 */
template <class M,
          class SGV,
          size_t s_r = 1,
          size_t s_rC = 1,
          class F = double,
          size_t r_r = s_r,
          size_t r_rC = s_rC,
          class RGV = SGV>
class LocalIntersectionOperatorFiniteDifferenceJacobianAssembler : public XT::Grid::IntersectionFunctor<SGV>
{
  static_assert(XT::LA::is_matrix<M>::value, "");
  static_assert(XT::Grid::is_view<SGV>::value, "");
  static_assert(XT::Grid::is_view<RGV>::value, "");
  static_assert(std::is_same<XT::Grid::extract_entity_t<SGV>, XT::Grid::extract_entity_t<RGV>>::value, "");

  using BaseType = XT::Grid::IntersectionFunctor<SGV>;
  using ThisType = LocalIntersectionOperatorFiniteDifferenceJacobianAssembler<M, SGV, s_r, s_rC, F, r_r, r_rC, RGV>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::I;
  using typename BaseType::IntersectionType;

  using MatrixType = M;
  using VectorType = XT::LA::vector_t<M>;
  using V = VectorType;

  using LocalIntersectionOperatorType = LocalIntersectionOperatorInterface<I, V, SGV, s_r, s_rC, F, r_r, r_rC, F, RGV>;
  using SourceSpaceType = typename LocalIntersectionOperatorType::DiscreteSourceType::SpaceType;
  using RangeSpaceType = typename LocalIntersectionOperatorType::LocalInsideRangeType::SpaceType;

  LocalIntersectionOperatorFiniteDifferenceJacobianAssembler(const SourceSpaceType& source_space,
                                                             const RangeSpaceType& range_space,
                                                             MatrixType& matrix,
                                                             const VectorType& source_vector,
                                                             const LocalIntersectionOperatorType& local_operator,
                                                             const XT::Common::Parameter& param = {},
                                                             const real_t<F> eps = 1e-7)
    : source_space_(source_space)
    , range_space_(range_space)
    , matrix_(matrix)
    , source_vector_(source_vector)
    , param_(param)
    , scaling_(param_.has_key("matrixoperator.scaling") ? param_.get("matrixoperator.scaling").at(0) : 1.)
    , eps_(eps)
    , source_(source_space_) // This is a full vector, and intended!
    , range_(range_space_) // This is a full vector, and intended!
    , local_source_inside_(source_.local_discrete_function())
    , local_source_outside_(source_.local_discrete_function())
    , local_range_inside_(range_.local_discrete_function())
    , local_range_outside_(range_.local_discrete_function())
    , local_op_(local_operator.with_source(source_))
  {
    source_.dofs().vector() = source_vector_;
  }

  LocalIntersectionOperatorFiniteDifferenceJacobianAssembler(const ThisType& other)
    : BaseType(other)
    , source_space_(other.source_space_)
    , range_space_(other.range_space_)
    , matrix_(other.matrix_)
    , source_vector_(other.source_vector_)
    , param_(other.param_)
    , scaling_(other.scaling_)
    , eps_(other.eps_)
    , source_(source_space_) // This is a full vector, and intended!
    , range_(range_space_) // This is a full vector, and intended!
    , local_source_inside_(source_.local_discrete_function())
    , local_source_outside_(source_.local_discrete_function())
    , local_range_inside_(range_.local_discrete_function())
    , local_range_outside_(range_.local_discrete_function())
    , local_op_(other.local_op_->with_source(source_))
  {
    source_.dofs().vector() = source_vector_;
  }

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    const bool treat_outside = intersection.neighbor();
    // some preparations
    local_source_inside_->bind(inside_element);
    local_range_inside_->bind(inside_element);
    source_space_.mapper().global_indices(inside_element, global_source_indices_inside_);
    range_space_.mapper().global_indices(inside_element, global_range_indices_inside_);
    const size_t local_source_inside_size = source_space_.mapper().local_size(inside_element);
    const size_t local_source_outside_size = treat_outside ? source_space_.mapper().local_size(outside_element) : 0;
    const size_t local_range_inside_size = range_space_.mapper().local_size(inside_element);
    const size_t local_range_outside_size = treat_outside ? range_space_.mapper().local_size(outside_element) : 0;
    if (range_DoFs_inside_.size() < local_range_inside_size)
      range_DoFs_inside_.resize(local_range_inside_size, 0);
    local_range_inside_->dofs().set_all(0);
    if (treat_outside) {
      local_source_outside_->bind(outside_element);
      local_range_outside_->bind(outside_element);
      source_space_.mapper().global_indices(outside_element, global_source_indices_outside_);
      range_space_.mapper().global_indices(outside_element, global_range_indices_outside_);
      if (range_DoFs_outside_.size() < local_range_outside_size)
        range_DoFs_outside_.resize(local_range_outside_size, 0);
      local_range_outside_->dofs().set_all(0);
    }
    // apply op as is, keep the result, clear local range
    local_op_->apply(intersection, *local_range_inside_, *local_range_outside_, param_);
    for (size_t ii = 0; ii < local_range_inside_size; ++ii)
      range_DoFs_inside_[ii] = local_range_inside_->dofs()[ii];
    local_range_inside_->dofs().set_all(0);
    if (treat_outside) {
      for (size_t ii = 0; ii < local_range_outside_size; ++ii)
        range_DoFs_outside_[ii] = local_range_outside_->dofs()[ii];
      local_range_outside_->dofs().set_all(0);
    }
    // loop over all inside source DoFs
    for (size_t jj = 0; jj < local_source_inside_size; ++jj) {
      // perturb source DoF
      const auto jjth_source_DoF = local_source_inside_->dofs()[jj];
      const auto eps = eps_ * (1. + std::abs(jjth_source_DoF));
      local_source_inside_->dofs()[jj] += eps;
      // apply op with perturbed source DoF
      local_op_->apply(intersection, *local_range_inside_, *local_range_outside_, param_);
      // observe perturbation in inside range DoFs
      for (size_t ii = 0; ii < local_range_inside_size; ++ii) {
        auto derivative = (local_range_inside_->dofs()[ii] - range_DoFs_inside_[ii]) / eps;
        if (XT::Common::FloatCmp::eq(derivative, eps))
          derivative = 0;
        matrix_.add_to_entry(
            global_range_indices_inside_[ii], global_source_indices_inside_[jj], scaling_ * derivative);
      }
      // observe perturbation in outside range DoFs
      if (treat_outside) {
        for (size_t ii = 0; ii < local_range_outside_size; ++ii) {
          auto derivative = (local_range_outside_->dofs()[ii] - range_DoFs_outside_[ii]) / eps;
          if (XT::Common::FloatCmp::eq(derivative, eps))
            derivative = 0;
          matrix_.add_to_entry(
              global_range_indices_outside_[ii], global_source_indices_inside_[jj], scaling_ * derivative);
        }
      }
      // restore source
      local_source_inside_->dofs()[jj] = jjth_source_DoF;
    }
    if (treat_outside) {
      // clear local range
      local_range_inside_->dofs().set_all(0);
      local_range_outside_->dofs().set_all(0);
      // loop over all outside source DoFs
      for (size_t jj = 0; jj < local_source_outside_size; ++jj) {
        // perturb source DoF
        const auto jjth_source_DoF = local_source_outside_->dofs()[jj];
        const auto eps = eps_ * (1. + std::abs(jjth_source_DoF));
        local_source_outside_->dofs()[jj] += eps;
        // apply op with perturbed source DoF
        local_op_->apply(intersection, *local_range_inside_, *local_range_outside_, param_);
        // observe perturbation in inside range DoFs
        for (size_t ii = 0; ii < local_range_inside_size; ++ii) {
          auto derivative = (local_range_inside_->dofs()[ii] - range_DoFs_inside_[ii]) / eps;
          if (XT::Common::FloatCmp::eq(derivative, eps))
            derivative = 0;
          matrix_.add_to_entry(
              global_range_indices_inside_[ii], global_source_indices_outside_[jj], scaling_ * derivative);
        }
        // observe perturbation in outside range DoFs
        for (size_t ii = 0; ii < local_range_outside_size; ++ii) {
          auto derivative = (local_range_outside_->dofs()[ii] - range_DoFs_outside_[ii]) / eps;
          if (XT::Common::FloatCmp::eq(derivative, eps))
            derivative = 0;
          matrix_.add_to_entry(
              global_range_indices_outside_[ii], global_source_indices_outside_[jj], scaling_ * derivative);
        }
        // restore source
        local_source_outside_->dofs()[jj] = jjth_source_DoF;
      }
    }
  } // ... apply_local(...)

private:
  const SourceSpaceType& source_space_;
  const RangeSpaceType& range_space_;
  MatrixType& matrix_;
  const VectorType& source_vector_;
  const XT::Common::Parameter param_;
  const double scaling_;
  const real_t<F> eps_;
  DiscreteFunction<V, SGV, s_r, s_rC, F> source_;
  DiscreteFunction<V, RGV, r_r, r_rC, F> range_;
  std::unique_ptr<LocalDiscreteFunction<V, SGV, s_r, s_rC, F>> local_source_inside_;
  std::unique_ptr<LocalDiscreteFunction<V, SGV, s_r, s_rC, F>> local_source_outside_;
  std::unique_ptr<LocalDiscreteFunction<V, RGV, r_r, r_rC, F>> local_range_inside_;
  std::unique_ptr<LocalDiscreteFunction<V, RGV, r_r, r_rC, F>> local_range_outside_;
  DynamicVector<size_t> global_source_indices_inside_;
  DynamicVector<size_t> global_source_indices_outside_;
  DynamicVector<size_t> global_range_indices_inside_;
  DynamicVector<size_t> global_range_indices_outside_;
  DynamicVector<F> range_DoFs_inside_;
  DynamicVector<F> range_DoFs_outside_;
  const std::unique_ptr<LocalIntersectionOperatorType> local_op_;
}; // class LocalIntersectionOperatorFiniteDifferenceJacobianAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_OPERATOR_FD_JACOBIAN_ASSEMBLERS_HH
