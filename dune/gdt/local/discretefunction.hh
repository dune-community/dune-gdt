// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2011 - 2018)
//   René Fritze     (2014, 2016, 2018)
//   René Milk       (2017)
//   Tobias Leibner  (2014, 2016 - 2017)

#ifndef DUNE_GDT_LOCAL_DISCRETEFUNCTION_HH
#define DUNE_GDT_LOCAL_DISCRETEFUNCTION_HH

#include <dune/xt/la/container/vector-interface.hh>
#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/grid-function.hh>

#include <dune/gdt/discretefunction/dof-vector.hh>
#include <dune/gdt/exceptions.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class Vector, class GridView, size_t range_dim = 1, size_t range_dim_cols = 1, class RangeField = double>
class ConstLocalDiscreteFunction
  : public XT::Functions::
        ElementFunctionInterface<XT::Grid::extract_entity_t<GridView>, range_dim, range_dim_cols, RangeField>
{
  // No need to check the rest, is done in SpaceInterface.
  static_assert(XT::LA::is_vector<Vector>::value, "");
  static_assert(range_dim_cols == 1, "The matrix-valued case requires updates to evaluate/jacobian/derivative!");

  using ThisType = ConstLocalDiscreteFunction;
  using BaseType = XT::Functions::
      ElementFunctionInterface<XT::Grid::extract_entity_t<GridView>, range_dim, range_dim_cols, RangeField>;

public:
  using SpaceType = SpaceInterface<GridView, range_dim, range_dim_cols, RangeField>;
  using ConstDofVectorType = ConstDofVector<Vector, GridView>;
  using ConstLocalDofVectorType = typename ConstDofVectorType::ConstLocalDofVectorType;
  using LocalBasisType = typename SpaceType::GlobalBasisType::LocalizedType;

  using BaseType::d;
  using BaseType::r;
  using BaseType::rC;
  using typename BaseType::DerivativeRangeReturnType;
  using typename BaseType::DerivativeRangeSelector;
  using typename BaseType::DerivativeRangeType;
  using typename BaseType::DomainType;
  using typename BaseType::DynamicDerivativeRangeType;
  using typename BaseType::DynamicRangeType;
  using typename BaseType::ElementType;
  using typename BaseType::R;
  using typename BaseType::RangeReturnType;
  using typename BaseType::RangeSelector;
  using typename BaseType::RangeType;
  using typename BaseType::SingleDerivativeRangeReturnType;
  using typename BaseType::SingleDerivativeRangeType;

  ConstLocalDiscreteFunction(const SpaceType& spc, const ConstDofVectorType& dof_vector)
    : BaseType()
    , space_(spc.copy())
    , space_is_fv_(space_->type() == GDT::SpaceType::finite_volume)
    , dof_vector_(dof_vector.localize())
    , basis_(space_->basis().localize())
    , basis_values_(space_is_fv_ ? 0 : space_->mapper().max_local_size())
    , dynamic_basis_values_(basis_values_.size())
    , basis_derivatives_(basis_values_.size())
    , dynamic_basis_derivatives_(basis_values_.size())
  {}

  virtual ~ConstLocalDiscreteFunction() = default;

protected:
  void post_bind(const ElementType& ent) override
  {
    basis_->bind(ent);
    basis_size_ = basis_->size();
    dof_vector_.bind(ent);
  }

public:
  int order(const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return basis_->order();
  }

  const SpaceType& space() const
  {
    return *space_;
  }

  const LocalBasisType& basis() const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return *basis_;
  }

  const ConstLocalDofVectorType& dofs() const
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return dof_vector_;
  }

  using BaseType::derivative;
  using BaseType::evaluate;
  using BaseType::jacobian;

  /**
   * \name ``These methods are required by XT::Functions::GridFunctionInterface.''
   * \{
   **/

  RangeReturnType evaluate(const DomainType& point_in_reference_element,
                           const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    RangeReturnType result(0);
    if (space_is_fv_) {
      for (size_t ii = 0; ii < r; ++ii)
        result[ii] = dof_vector_[ii];
    } else {
      basis_->evaluate(point_in_reference_element, basis_values_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result.axpy(dof_vector_[ii], basis_values_[ii]);
    }
    return result;
  } // ... evaluate(...)

  DerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                     const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    DerivativeRangeReturnType result(0);
    if (space_is_fv_) {
      return result;
    } else {
      basis_->jacobians(point_in_reference_element, basis_derivatives_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result.axpy(dof_vector_[ii], basis_derivatives_[ii]);
    }
    return result;
  } // ... jacobian(...)

  DerivativeRangeReturnType derivative(const std::array<size_t, d>& alpha,
                                       const DomainType& point_in_reference_element,
                                       const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    DUNE_THROW_IF(space_is_fv_,
                  Exceptions::discrete_function_error,
                  "arbitrary derivatives are not supported by the local finite elements!\n\n"
                      << "alpha = " << alpha << "\n"
                      << "point_in_reference_element = " << point_in_reference_element);
    DerivativeRangeReturnType result(0);
    for (size_t jj = 0; jj < d; ++jj) {
      if (alpha[jj] == 0) {
        for (size_t ii = 0; ii < r; ++ii)
          result[ii][jj] = dof_vector_[ii];
      } else
        DUNE_THROW(Exceptions::discrete_function_error,
                   "arbitrary derivatives are not supported by the local finite elements!\n\n"
                       << "alpha = " << alpha << "\n"
                       << "point_in_reference_element = " << point_in_reference_element);
    }
    return result;
  } // ... derivative(...)

  /**
   * \}
   * \name ``These methods are default implemented in XT::Functions::GridFunctionInterface and are overridden
   *         for improved performance.''
   * \{
   **/

  void evaluate(const DomainType& point_in_reference_element,
                DynamicRangeType& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    RangeSelector::ensure_size(result);
    if (space_is_fv_) {
      for (size_t ii = 0; ii < r; ++ii)
        result[ii] = dof_vector_[ii];
    } else {
      result *= 0;
      basis_->evaluate(point_in_reference_element, dynamic_basis_values_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result.axpy(dof_vector_[ii], dynamic_basis_values_[ii]);
    }
  } // ... evaluate(...)

  void jacobian(const DomainType& point_in_reference_element,
                DynamicDerivativeRangeType& result,
                const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    DerivativeRangeSelector::ensure_size(result);
    result *= 0;
    if (space_is_fv_) {
      return;
    } else {
      basis_->jacobians(point_in_reference_element, dynamic_basis_derivatives_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result.axpy(dof_vector_[ii], dynamic_basis_derivatives_[ii]);
    }
  } // ... jacobian(...)

  void derivative(const std::array<size_t, d>& alpha,
                  const DomainType& point_in_reference_element,
                  DynamicDerivativeRangeType& result,
                  const XT::Common::Parameter& /*param*/ = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    DUNE_THROW_IF(space_is_fv_,
                  Exceptions::discrete_function_error,
                  "arbitrary derivatives are not supported by the local finite elements!\n\n"
                      << "alpha = " << alpha << "\n"
                      << "point_in_reference_element = " << point_in_reference_element);
    DerivativeRangeSelector::ensure_size(result);
    result *= 0;
    for (size_t jj = 0; jj < d; ++jj) {
      if (alpha[jj] == 0) {
        for (size_t ii = 0; ii < r; ++ii)
          result[ii][jj] = dof_vector_[ii];
      } else
        DUNE_THROW(Exceptions::discrete_function_error,
                   "arbitrary derivatives are not supported by the local finite elements!\n\n"
                       << "alpha = " << alpha << "\n"
                       << "point_in_reference_element = " << point_in_reference_element);
    }
  } // ... derivative(...)

  /**
   * \}
   * \name ``These methods (used to access an individual range dimension) are default implemented in
   *         XT::Functions::GridFunctionInterface and are implemented for improved performance.''
   * \{
   **/

  R evaluate(const DomainType& point_in_reference_element,
             const size_t row,
             const size_t col = 0,
             const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    this->assert_correct_dims(row, col, "evaluate");
    if (space_is_fv_) {
      return dof_vector_[row];
    } else {
      R result(0);
      basis_->evaluate(point_in_reference_element, basis_values_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result += dof_vector_[ii] * basis_values_[ii][row];
      return result;
    }
  } // ... evaluate(...)

  SingleDerivativeRangeReturnType jacobian(const DomainType& point_in_reference_element,
                                           const size_t row,
                                           const size_t col = 0,
                                           const XT::Common::Parameter& param = {}) const override final
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    this->assert_correct_dims(row, col, "jacobian");
    if (space_is_fv_) {
      return 0;
    } else {
      SingleDerivativeRangeReturnType result(0);
      basis_->jacobians(point_in_reference_element, basis_derivatives_, param);
      for (size_t ii = 0; ii < basis_size_; ++ii)
        result.axpy(dof_vector_[ii], basis_derivatives_[ii][row]);
      return result;
    }
  } // ... jacobian(...)

  /// \}

private:
  std::unique_ptr<const SpaceType> space_;
  const bool space_is_fv_;
  ConstLocalDofVectorType dof_vector_;
  std::unique_ptr<LocalBasisType> basis_;
  mutable size_t basis_size_;
  mutable std::vector<RangeType> basis_values_;
  mutable std::vector<DynamicRangeType> dynamic_basis_values_;
  mutable std::vector<DerivativeRangeType> basis_derivatives_;
  mutable std::vector<DynamicDerivativeRangeType> dynamic_basis_derivatives_;
}; // class ConstLocalDiscreteFunction


template <class V, class GV, size_t r = 1, size_t rC = 1, class R = double>
class LocalDiscreteFunction : public ConstLocalDiscreteFunction<V, GV, r, rC, R>
{
  using ThisType = LocalDiscreteFunction;
  using BaseType = ConstLocalDiscreteFunction<V, GV, r, rC, R>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::SpaceType;

  using DofVectorType = DofVector<V, GV>;
  using LocalDofVectorType = typename DofVectorType::LocalDofVectorType;

  LocalDiscreteFunction(const SpaceType& spc, DofVectorType& dof_vector)
    : BaseType(spc, dof_vector)
    , dof_vector_(dof_vector.localize())
  {}

  LocalDiscreteFunction(const SpaceType& spc, DofVectorType& dof_vector, const ElementType& ent)
    : BaseType(spc, dof_vector, ent)
    , dof_vector_(dof_vector.localize(ent))
  {}

protected:
  void post_bind(const ElementType& ent) override final
  {
    BaseType::post_bind(ent);
    dof_vector_.bind(ent);
  }

public:
  LocalDofVectorType& dofs()
  {
    DUNE_THROW_IF(!this->is_bound_, Exceptions::not_bound_to_an_element_yet, "");
    return dof_vector_;
  }

private:
  LocalDofVectorType dof_vector_;
}; // class LocalDiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_DISCRETEFUNCTION_HH
