// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016 - 2017)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_LOCAL_ASSEMLBER_HH
#define DUNE_GDT_LOCAL_ASSEMLBER_HH

#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/grid/walker/apply-on.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/la/container/interfaces.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/local/discretefunction.hh>
#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class LocalVolumeTwoFormType>
class LocalVolumeTwoFormAssembler
{
  static_assert(is_local_volume_twoform<LocalVolumeTwoFormType>::value,
                "LocalVolumeTwoFormType has to be derived from LocalVolumeTwoFormInterface!");

public:
  explicit LocalVolumeTwoFormAssembler(const LocalVolumeTwoFormType& local_twoform)
    : local_volume_twoform_(local_twoform)
  {
  }

  /**
   *  \tparam T           Traits of the SpaceInterface implementation, representing the type of test_space
   *  \tparam A           Traits of the SpaceInterface implementation, representing the type of ansatz_space
   *  \tparam *d          dimDomain of test_space (* == T) or ansatz_space (* == A)
   *  \tparam *r          dimRange of test_space (* == T) or ansatz_space (* == A)
   *  \tparam *rC         dimRangeCols of test_space (* == T) or ansatz_space (* == A)
   *  \tparam EntityType  A model of Dune::Entity< 0 >
   *  \tparam M           Traits of the Dune::XT::LA::Container::MatrixInterface implementation, representing the
   * type of global_matrix
   *  \tparam R           RangeFieldType, i.e. double
   */
  template <class T,
            size_t Td,
            size_t Tr,
            size_t TrC,
            class A,
            size_t Ad,
            size_t Ar,
            size_t ArC,
            class EntityType,
            class M,
            class R>
  void assemble(const SpaceInterface<T, Td, Tr, TrC>& test_space,
                const SpaceInterface<A, Ad, Ar, ArC>& ansatz_space,
                const EntityType& entity,
                XT::LA::MatrixInterface<M, R>& global_matrix) const
  {
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    Dune::DynamicMatrix<R> local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_volume_twoform_.apply2(test_base, ansatz_base, local_matrix);
    // write local matrix to global
    const auto global_row_indices =
        test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    const auto global_col_indices =
        ansatz_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_row_indices.size() == rows);
    assert(global_col_indices.size() == cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& local_matrix_row = local_matrix[ii];
      const size_t global_ii = global_row_indices[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t global_jj = global_col_indices[jj];
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix_row[jj]);
      }
    } // write local matrix to global
  } // ... assemble(...)

private:
  const LocalVolumeTwoFormType& local_volume_twoform_;
}; // class LocalVolumeTwoFormAssembler


template <class GridViewImp,
          class LocalVolumeTwoFormType,
          class TestFunctionType,
          class AnsatzFunctionType,
          class FieldType>
class LocalVolumeTwoFormAccumulator : public XT::Grid::internal::Codim0ReturnObject<GridViewImp, FieldType>
{
  static_assert(std::is_base_of<LocalVolumeTwoFormInterface<typename LocalVolumeTwoFormType::Traits>,
                                LocalVolumeTwoFormType>::value,
                "LocalVolumeTwoFormType has to be derived from LocalVolumeTwoFormInterface!");
  static_assert(XT::Functions::is_localizable_function<TestFunctionType>::value,
                "TestFunctionType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(XT::Functions::is_localizable_function<AnsatzFunctionType>::value,
                "AnsatzFunctionType has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef LocalVolumeTwoFormAccumulator<GridViewImp,
                                        LocalVolumeTwoFormType,
                                        TestFunctionType,
                                        AnsatzFunctionType,
                                        FieldType>
      ThisType;
  typedef XT::Grid::internal::Codim0ReturnObject<GridViewImp, FieldType> BaseType;

public:
  typedef typename BaseType::GridViewType GridViewType;
  typedef typename BaseType::EntityType EntityType;

  LocalVolumeTwoFormAccumulator(const GridViewType& grd_vw,
                                const LocalVolumeTwoFormType& local_op,
                                const TestFunctionType& test_function,
                                const AnsatzFunctionType& ansatz_function,
                                const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where)
    : grid_view_(grd_vw)
    , local_operator_(local_op)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , result_(0)
    , finalized_(false)
    , where_(where)
  {
  }

  LocalVolumeTwoFormAccumulator(const ThisType& other) = default;
  virtual ~LocalVolumeTwoFormAccumulator() = default;

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const override final
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual FieldType compute_locally(const EntityType& entity) override final
  {
    DynamicMatrix<FieldType> local_twoform_result(1, 1, 0.); // \todo: make mutable member, after SMP refactor
    this->local_operator_.apply2(
        *test_function_.local_function(entity), *ansatz_function_.local_function(entity), local_twoform_result);
    return local_twoform_result[0][0];
  } // ... compute_locally(...)

  virtual void apply_local(const EntityType& entity) override
  {
    *result_ += compute_locally(entity);
  }

  virtual void finalize() override
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  virtual FieldType result() const override final
  {
    if (!finalized_)
      DUNE_THROW(XT::Common::Exceptions::you_are_using_this_wrong, "Call finalize() first!");
    return finalized_result_;
  }

private:
  const GridViewType& grid_view_;
  const LocalVolumeTwoFormType& local_operator_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  Dune::XT::Common::PerThreadValue<FieldType> result_;
  bool finalized_;
  const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where_;
  FieldType finalized_result_;
}; // class LocalVolumeTwoFormAccumulator


template <class GridViewType, class LocalOperatorType, class SourceType, class RangeType>
class LocalOperatorApplicator : public XT::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(is_local_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim0Object<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  LocalOperatorApplicator(const GridViewType& grid_view,
                          const LocalOperatorType& local_operator,
                          const SourceType& source,
                          RangeType& range,
                          const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.apply(source_, *range_.local_discrete_function(entity));
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where_;
}; // class LocalOperatorApplicator


template <class GridViewType, class LocalOperatorType, class SourceType, class RangeType>
class LocalOperatorJacobianAssembler : public XT::Grid::internal::Codim0Object<GridViewType>
{
  static_assert(is_local_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim0Object<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;

  LocalOperatorJacobianAssembler(const GridViewType& grid_view,
                                 const LocalOperatorType& local_operator,
                                 const SourceType& x,
                                 const SourceType& source,
                                 RangeType& range,
                                 const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , x_(x)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const EntityType& entity) const
  {
    return where_.apply_on(grid_view, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.assemble_jacobian(x_, source_, *range_.local_discrete_function(entity));
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& x_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichEntity<GridViewType>& where_;
}; // class LocalOperatorJacobianAssembler


template <class LocalCouplingTwoFormType>
class LocalCouplingTwoFormAssembler
{
  static_assert(is_local_coupling_twoform<LocalCouplingTwoFormType>::value,
                "LocalCouplingTwoFormType has to be derived from LocalCouplingTwoFormInterface!");

public:
  explicit LocalCouplingTwoFormAssembler(const LocalCouplingTwoFormType& local_twoform)
    : local_coupling_twoform_(local_twoform)
  {
  }

  template <class TE,
            size_t TEd,
            size_t TEr,
            size_t TErC,
            class AE,
            size_t AEd,
            size_t AEr,
            size_t AErC,
            class TN,
            size_t TNd,
            size_t TNr,
            size_t TNrC,
            class AN,
            size_t ANd,
            size_t ANr,
            size_t ANrC,
            class IntersectionType,
            class MEE,
            class MNN,
            class MEN,
            class MNE,
            class R>
  void assemble(const SpaceInterface<TE, TEd, TEr, TErC>& test_space_en,
                const SpaceInterface<AE, AEd, AEr, AErC>& ansatz_space_en,
                const SpaceInterface<TN, TNd, TNr, TNrC>& test_space_ne,
                const SpaceInterface<AN, ANd, ANr, ANrC>& ansatz_space_ne,
                const IntersectionType& intersection,
                XT::LA::MatrixInterface<MEE, R>& global_matrix_en_en,
                XT::LA::MatrixInterface<MNN, R>& global_matrix_ne_ne,
                XT::LA::MatrixInterface<MEN, R>& global_matrix_en_ne,
                XT::LA::MatrixInterface<MNE, R>& global_matrix_ne_en) const
  {
    assert(global_matrix_en_en.rows() >= test_space_en.mapper().size());
    assert(global_matrix_en_en.cols() >= ansatz_space_en.mapper().size());
    assert(global_matrix_ne_ne.rows() >= test_space_ne.mapper().size());
    assert(global_matrix_en_en.cols() >= ansatz_space_ne.mapper().size());
    assert(global_matrix_en_ne.rows() >= test_space_en.mapper().size());
    assert(global_matrix_en_ne.cols() >= ansatz_space_ne.mapper().size());
    assert(global_matrix_ne_en.rows() >= test_space_ne.mapper().size());
    assert(global_matrix_ne_en.cols() >= ansatz_space_en.mapper().size());
    const auto entity = intersection.inside();
    const auto neighbor = intersection.outside();
    // prepare
    const size_t rows_en = test_space_en.mapper().numDofs(entity);
    const size_t cols_en = ansatz_space_en.mapper().numDofs(entity);
    const size_t rows_ne = test_space_ne.mapper().numDofs(neighbor);
    const size_t cols_ne = ansatz_space_ne.mapper().numDofs(neighbor);
    Dune::DynamicMatrix<R> local_matrix_en_en(rows_en, cols_en, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> local_matrix_ne_ne(rows_ne, cols_ne, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> local_matrix_en_ne(rows_en, cols_ne, 0.); // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<R> local_matrix_ne_en(rows_ne, cols_en, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base_en = test_space_en.base_function_set(entity);
    const auto ansatz_base_en = ansatz_space_en.base_function_set(entity);
    const auto test_base_ne = test_space_ne.base_function_set(neighbor);
    const auto ansatz_base_ne = ansatz_space_ne.base_function_set(neighbor);
    local_coupling_twoform_.apply2(test_base_en,
                                   ansatz_base_en,
                                   test_base_ne,
                                   ansatz_base_ne,
                                   intersection,
                                   local_matrix_en_en,
                                   local_matrix_ne_ne,
                                   local_matrix_en_ne,
                                   local_matrix_ne_en);
    // write local matrix to global
    const auto global_row_indices_en =
        test_space_en.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    const auto global_col_indices_en =
        ansatz_space_en.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    const auto global_row_indices_ne =
        test_space_ne.mapper().globalIndices(neighbor); // \todo: make mutable member, after SMP refactor
    const auto global_col_indices_ne =
        ansatz_space_ne.mapper().globalIndices(neighbor); // \todo: make mutable member, after SMP refactor
    assert(global_row_indices_en.size() == rows_en);
    assert(global_col_indices_en.size() == cols_en);
    assert(global_row_indices_ne.size() == rows_ne);
    assert(global_col_indices_ne.size() == cols_ne);
    for (size_t ii = 0; ii < rows_en; ++ii) {
      const auto& local_matrix_en_en_row = local_matrix_en_en[ii];
      const auto& local_matrix_en_ne_row = local_matrix_en_ne[ii];
      const size_t global_ii = global_row_indices_en[ii];
      for (size_t jj = 0; jj < cols_en; ++jj) {
        const size_t global_jj = global_col_indices_en[jj];
        global_matrix_en_en.add_to_entry(global_ii, global_jj, local_matrix_en_en_row[jj]);
      }
      for (size_t jj = 0; jj < cols_ne; ++jj) {
        const size_t global_jj = global_col_indices_ne[jj];
        global_matrix_en_ne.add_to_entry(global_ii, global_jj, local_matrix_en_ne_row[jj]);
      }
    }
    for (size_t ii = 0; ii < rows_ne; ++ii) {
      const auto& local_matrix_ne_en_row = local_matrix_ne_en[ii];
      const auto& local_matrix_ne_ne_row = local_matrix_ne_ne[ii];
      const size_t global_ii = global_row_indices_ne[ii];
      for (size_t jj = 0; jj < cols_en; ++jj) {
        const size_t global_jj = global_col_indices_en[jj];
        global_matrix_ne_en.add_to_entry(global_ii, global_jj, local_matrix_ne_en_row[jj]);
      }
      for (size_t jj = 0; jj < cols_ne; ++jj) {
        const size_t global_jj = global_col_indices_ne[jj];
        global_matrix_ne_ne.add_to_entry(global_ii, global_jj, local_matrix_ne_ne_row[jj]);
      }
    }
  } // ... assemble(...)

  template <class TE,
            size_t TEd,
            size_t TEr,
            size_t TErC,
            class AE,
            size_t AEd,
            size_t AEr,
            size_t AErC,
            class TN,
            size_t TNd,
            size_t TNr,
            size_t TNrC,
            class AN,
            size_t ANd,
            size_t ANr,
            size_t ANrC,
            class IntersectionType,
            class M,
            class R>
  void assemble(const SpaceInterface<TE, TEd, TEr, TErC>& test_space_en,
                const SpaceInterface<AE, AEd, AEr, AErC>& ansatz_space_en,
                const SpaceInterface<TN, TNd, TNr, TNrC>& test_space_ne,
                const SpaceInterface<AN, ANd, ANr, ANrC>& ansatz_space_ne,
                const IntersectionType& intersection,
                XT::LA::MatrixInterface<M, R>& global_matrix) const
  {
    assemble(test_space_en,
             ansatz_space_en,
             test_space_ne,
             ansatz_space_ne,
             intersection,
             global_matrix,
             global_matrix,
             global_matrix,
             global_matrix);
  }

private:
  const LocalCouplingTwoFormType& local_coupling_twoform_;
}; // class LocalCouplingTwoFormAssembler


template <class GridViewType, class LocalOperatorType, class SourceType, class RangeType>
class LocalCouplingOperatorApplicator : public XT::Grid::internal::Codim1Object<GridViewType>
{
  static_assert(is_local_coupling_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalCouplingOperatorApplicator(const GridViewType& grid_view,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridViewType>& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_view, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity,
                           const double t = 0)
  {
    local_operator_.apply(source_,
                          intersection,
                          *range_.local_discrete_function(inside_entity),
                          *range_.local_discrete_function(outside_entity),
                          t);
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridViewType>& where_;
}; // class LocalCouplingOperatorApplicator


template <class LocalBoundaryTwoFormType>
class LocalBoundaryTwoFormAssembler
{
  static_assert(is_local_boundary_twoform<LocalBoundaryTwoFormType>::value, "");

public:
  explicit LocalBoundaryTwoFormAssembler(const LocalBoundaryTwoFormType& local_twoform)
    : local_twoform_(local_twoform)
  {
  }

  template <class T,
            size_t Td,
            size_t Tr,
            size_t TrC,
            class A,
            size_t Ad,
            size_t Ar,
            size_t ArC,
            class IntersectionType,
            class M,
            class R>
  void assemble(const SpaceInterface<T, Td, Tr, TrC>& test_space,
                const SpaceInterface<A, Ad, Ar, ArC>& ansatz_space,
                const IntersectionType& intersection,
                XT::LA::MatrixInterface<M, R>& global_matrix) const
  {
    const auto entity = intersection.inside();
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    Dune::DynamicMatrix<R> local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_twoform_.apply2(test_base, ansatz_base, intersection, local_matrix);
    // write local matrix to global
    const auto global_row_indices =
        test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    const auto global_col_indices =
        ansatz_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_row_indices.size() == rows);
    assert(global_col_indices.size() == cols);
    for (size_t ii = 0; ii < rows; ++ii) {
      const auto& local_matrix_row = local_matrix[ii];
      const size_t global_ii = global_row_indices[ii];
      for (size_t jj = 0; jj < cols; ++jj) {
        const size_t global_jj = global_col_indices[jj];
        global_matrix.add_to_entry(global_ii, global_jj, local_matrix_row[jj]);
      }
    } // write local matrix to global
  } // ... assemble(...)

private:
  const LocalBoundaryTwoFormType& local_twoform_;
}; // class LocalBoundaryTwoFormAssembler


template <class GridViewType, class LocalOperatorType, class SourceType, class RangeType>
class LocalBoundaryOperatorApplicator : public XT::Grid::internal::Codim1Object<GridViewType>
{
  static_assert(is_local_boundary_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridViewType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalBoundaryOperatorApplicator(const GridViewType& grid_view,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridViewType>& where)
    : grid_view_(grid_view)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridViewType& grid_view, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_view, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& /*outside_entity*/,
                           const double t = 0)
  {
    local_operator_.apply(source_, intersection, *range_.local_discrete_function(inside_entity), t);
  }

private:
  const GridViewType& grid_view_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridViewType>& where_;
}; // class LocalBoundaryOperatorApplicator


template <class LocalVolumeFunctionalType>
class LocalVolumeFunctionalAssembler
{
  static_assert(is_local_volume_functional<LocalVolumeFunctionalType>::value,
                "LocalVolumeFunctionalType has to be derived from LocalVolumeFunctionalInterface!");

public:
  explicit LocalVolumeFunctionalAssembler(const LocalVolumeFunctionalType& local_volume_functional)
    : local_volume_functional_(local_volume_functional)
  {
  }

  /**
   *  \tparam S          Traits of the SpaceInterface implementation, representing the type of test_space
   *  \tparam d          dimDomain of test_space
   *  \tparam r          dimRange of test_space
   *  \tparam rC         dimRangeCols of test_space
   *  \tparam EntityType A model of Dune::Entity< 0 >
   *  \tparam V          Traits of the Dune::XT::LA::Container::VectorInterface implementation, representing the type
   * of global_vector
   *  \tparam R          RangeFieldType, i.e. double
   */
  template <class S, size_t d, size_t r, size_t rC, class EntityType, class V, class R>
  void assemble(const SpaceInterface<S, d, r, rC>& test_space,
                const EntityType& entity,
                XT::LA::VectorInterface<V, R>& global_vector) const
  {
    // prepare
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector<R> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_volume_functional_.apply(test_basis, local_vector);
    // write local vector to global
    const auto global_indices =
        test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalVolumeFunctionalType& local_volume_functional_;
}; // class LocalVolumeFunctionalAssembler


template <class LocalFunctionalType>
class LocalFaceFunctionalAssembler
{
  static_assert(
      std::is_base_of<LocalFaceFunctionalInterface<typename LocalFunctionalType::Traits>, LocalFunctionalType>::value,
      "LocalFunctionalType has to be derived from LocalFaceFunctionalInterface!");

public:
  explicit LocalFaceFunctionalAssembler(const LocalFunctionalType& local_face_functional)
    : local_face_functional_(local_face_functional)
  {
  }

  template <class T, size_t d, size_t r, size_t rC, class IntersectionType, class V, class R>
  void assemble(const SpaceInterface<T, d, r, rC>& test_space,
                const IntersectionType& intersection,
                XT::LA::VectorInterface<V, R>& global_vector) const
  {
    // prepare
    const auto entity = intersection.inside();
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector<R> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_face_functional_.apply(test_basis, intersection, local_vector);
    // write local vector to global
    const auto global_indices =
        test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalFunctionalType& local_face_functional_;
}; // class LocalFaceFunctionalAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMLBER_HH
