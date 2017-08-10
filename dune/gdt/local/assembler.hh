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
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


/**
 * \todo \attention Rename LocalVolumeTwoFormAssemblerFunctor -> LocalVolumeTwoFormAssembler after removing this class!
 */
template <class TestSpace, class Matrix, class AnsatzSpace = TestSpace>
class DUNE_DEPRECATED_MSG("Use LocalVolumeTwoFormAssemblerFunctor instead (13.05.2017)!") LocalVolumeTwoFormAssembler
{
  static_assert(is_space<TestSpace>::value, "");
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(is_space<AnsatzSpace>::value, "");
  static_assert(std::is_same<typename TestSpace::EntityType, typename AnsatzSpace::EntityType>::value, "");

public:
  typedef TestSpace TestSpaceType;
  typedef AnsatzSpace AnsatzSpaceType;
  typedef typename TestSpaceType::EntityType EntityType;
  typedef Matrix MatrixType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalVolumeTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                      typename AnsatzSpaceType::BaseFunctionSetType,
                                      FieldType>
      LocalVolumeTwoFormType;

  explicit LocalVolumeTwoFormAssembler(const LocalVolumeTwoFormType& local_twoform)
    : local_volume_twoform_(local_twoform)
  {
  }

  void assemble(const TestSpaceType& test_space,
                const AnsatzSpaceType& ansatz_space,
                const EntityType& entity,
                MatrixType& global_matrix) const
  {
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    DynamicMatrix<FieldType> local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
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
            class R>
  void assemble_entitywise(const SpaceInterface<T, Td, Tr, TrC>& test_space,
                           const SpaceInterface<A, Ad, Ar, ArC>& ansatz_space,
                           const EntityType& entity,
                           std::vector<Dune::DynamicMatrix<R>>& global_matrix) const
  {
// prepare
#ifndef NDEBUG
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
#endif
    auto& local_matrix = global_matrix[test_space.grid_view().indexSet().index(entity)];
    assert(local_matrix.N() == rows && local_matrix.M() == cols);
    local_matrix *= 0.;

    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_volume_twoform_.apply2(test_base, ansatz_base, local_matrix);
  } // ... assemble_entitywise(...)

private:
  const LocalVolumeTwoFormType& local_volume_twoform_;
}; // class LocalVolumeTwoFormAssembler


/**
 * \todo \attention Rename LocalVolumeTwoFormAccumulatorFunctor -> LocalVolumeTwoFormAccumulator after removing this
 *                  class!
 */
template <class GridLayerImp,
          class TestFunctionType,
          class AnsatzFunctionType,
          class FieldType = typename TestFunctionType::RangeFieldType>
class DUNE_DEPRECATED_MSG("User LocalVolumeTwoFormAccumulatorFunctor instead (29.05.2017)!")
    LocalVolumeTwoFormAccumulator : public XT::Grid::internal::Codim0ReturnObject<GridLayerImp, FieldType>
{
  static_assert(XT::Functions::is_localizable_function<TestFunctionType>::value,
                "TestFunctionType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(XT::Functions::is_localizable_function<AnsatzFunctionType>::value,
                "AnsatzFunctionType has to be derived from XT::Functions::LocalizableFunctionInterface!");

  typedef XT::Grid::internal::Codim0ReturnObject<GridLayerImp, FieldType> BaseType;

public:
  typedef LocalVolumeTwoFormInterface<typename TestFunctionType::LocalfunctionType,
                                      typename AnsatzFunctionType::LocalfunctionType,
                                      FieldType>
      LocalVolumeTwoFormType;
  typedef typename BaseType::GridLayerType GridLayerType;
  typedef typename BaseType::EntityType EntityType;

  LocalVolumeTwoFormAccumulator(const GridLayerType& grd_layr,
                                const LocalVolumeTwoFormType& local_op,
                                const TestFunctionType& test_function,
                                const AnsatzFunctionType& ansatz_function,
                                const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where)
    : grid_layer_(grd_layr)
    , local_operator_(local_op)
    , test_function_(test_function)
    , ansatz_function_(ansatz_function)
    , result_(0)
    , finalized_(false)
    , where_(where)
  {
  }

// disable deprecation warning that occurs even if this class is not used
#include <dune/xt/common/disable_warnings.hh>
  LocalVolumeTwoFormAccumulator(const LocalVolumeTwoFormAccumulator& other) = default;
#include <dune/xt/common/reenable_warnings.hh>

  virtual ~LocalVolumeTwoFormAccumulator() = default;

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const override final
  {
    return where_.apply_on(grid_layer, entity);
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
      finalized_result_ = grid_layer_.comm().sum(finalized_result_);
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
  const GridLayerType& grid_layer_;
  const LocalVolumeTwoFormType& local_operator_;
  const TestFunctionType& test_function_;
  const AnsatzFunctionType& ansatz_function_;
  Dune::XT::Common::PerThreadValue<FieldType> result_;
  bool finalized_;
  const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where_;
  FieldType finalized_result_;
}; // class LocalVolumeTwoFormAccumulator


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalOperatorApplicator : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_local_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;

  LocalOperatorApplicator(const GridLayerType& grid_layer,
                          const LocalOperatorType& local_operator,
                          const SourceType& source,
                          RangeType& range,
                          const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where)
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const EntityType& entity) const
  {
    return where_.apply_on(grid_layer, entity);
  }

  virtual void apply_local(const EntityType& entity)
  {
    local_operator_.apply(source_, *range_.local_discrete_function(entity));
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichEntity<GridLayerType>& where_;
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


/**
 * \todo \attention Rename LocalCouplingTwoFormAssemblerFunctor -> LocalCouplingTwoFormAssembler after removing this
 *                  class!
 */
template <class TestSpace,
          class Intersection,
          class Matrix,
          class AnsatzSpace = TestSpace,
          class OuterTestSpace = TestSpace,
          class OuterAnsatzSpace = AnsatzSpace>
class DUNE_DEPRECATED_MSG("Use LocalCouplingTwoFormAssemblerFunctor instead (13.05.2017)!")
    LocalCouplingTwoFormAssembler
{
  static_assert(is_space<TestSpace>::value, "");
  static_assert(is_space<AnsatzSpace>::value, "");
  static_assert(is_space<OuterTestSpace>::value, "");
  static_assert(is_space<OuterAnsatzSpace>::value, "");
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

public:
  typedef TestSpace TestSpaceType;
  typedef AnsatzSpace AnsatzSpaceType;
  typedef OuterTestSpace OuterTestSpaceType;
  typedef OuterAnsatzSpace OuterAnsatzSpaceType;
  typedef Matrix MatrixType;
  typedef Intersection IntersectionType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalCouplingTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        IntersectionType,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        typename OuterTestSpaceType::BaseFunctionSetType,
                                        typename OuterAnsatzSpaceType::BaseFunctionSetType,
                                        FieldType>
      LocalCouplingTwoFormType;

  explicit LocalCouplingTwoFormAssembler(const LocalCouplingTwoFormType& local_twoform)
    : local_coupling_twoform_(local_twoform)
  {
  }

  void assemble(const TestSpaceType& inner_test_space,
                const AnsatzSpaceType& inner_ansatz_space,
                const OuterTestSpaceType& outer_test_space,
                const OuterAnsatzSpaceType& outer_ansatz_space,
                const IntersectionType& intersection,
                MatrixType& global_matrix_in_in,
                MatrixType& global_matrix_out_out,
                MatrixType& global_matrix_in_out,
                MatrixType& global_matrix_out_in) const
  {
    assert(global_matrix_in_in.rows() >= inner_test_space.mapper().size());
    assert(global_matrix_in_in.cols() >= inner_ansatz_space.mapper().size());
    assert(global_matrix_out_out.rows() >= outer_test_space.mapper().size());
    assert(global_matrix_out_out.cols() >= outer_ansatz_space.mapper().size());
    assert(global_matrix_in_out.rows() >= inner_test_space.mapper().size());
    assert(global_matrix_in_out.cols() >= outer_ansatz_space.mapper().size());
    assert(global_matrix_out_in.rows() >= outer_test_space.mapper().size());
    assert(global_matrix_out_in.cols() >= inner_ansatz_space.mapper().size());
    const auto entity = intersection.inside();
    const auto neighbor = intersection.outside();
    // prepare
    const size_t rows_in = inner_test_space.mapper().numDofs(entity);
    const size_t cols_in = inner_ansatz_space.mapper().numDofs(entity);
    const size_t rows_out = outer_test_space.mapper().numDofs(neighbor);
    const size_t cols_out = outer_ansatz_space.mapper().numDofs(neighbor);
    // \todo: make mutable member, after SMP refactor
    DynamicMatrix<FieldType> local_matrix_in_in(rows_in, cols_in, 0.);
    DynamicMatrix<FieldType> local_matrix_out_out(rows_out, cols_out, 0.);
    DynamicMatrix<FieldType> local_matrix_in_out(rows_in, cols_out, 0.);
    DynamicMatrix<FieldType> local_matrix_out_in(rows_out, cols_in, 0.);
    // apply local two-form
    const auto test_base_in = inner_test_space.base_function_set(entity);
    const auto ansatz_base_in = inner_ansatz_space.base_function_set(entity);
    const auto test_base_out = outer_test_space.base_function_set(neighbor);
    const auto ansatz_base_out = outer_ansatz_space.base_function_set(neighbor);
    local_coupling_twoform_.apply2(test_base_in,
                                   ansatz_base_in,
                                   test_base_out,
                                   ansatz_base_out,
                                   intersection,
                                   local_matrix_in_in,
                                   local_matrix_out_out,
                                   local_matrix_in_out,
                                   local_matrix_out_in);
    // write local matrix to global
    // \todo: make mutable member, after SMP refactor
    const auto global_row_indices_in = inner_test_space.mapper().globalIndices(entity);
    const auto global_col_indices_in = inner_ansatz_space.mapper().globalIndices(entity);
    const auto global_row_indices_out = outer_test_space.mapper().globalIndices(neighbor);
    const auto global_col_indices_out = outer_ansatz_space.mapper().globalIndices(neighbor);
    assert(global_row_indices_in.size() == rows_in);
    assert(global_col_indices_in.size() == cols_in);
    assert(global_row_indices_out.size() == rows_out);
    assert(global_col_indices_out.size() == cols_out);
    for (size_t ii = 0; ii < rows_in; ++ii) {
      const auto& local_matrix_in_in_row = local_matrix_in_in[ii];
      const auto& local_matrix_in_out_row = local_matrix_in_out[ii];
      const size_t global_ii = global_row_indices_in[ii];
      for (size_t jj = 0; jj < cols_in; ++jj) {
        const size_t global_jj = global_col_indices_in[jj];
        global_matrix_in_in.add_to_entry(global_ii, global_jj, local_matrix_in_in_row[jj]);
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        const size_t global_jj = global_col_indices_out[jj];
        global_matrix_in_out.add_to_entry(global_ii, global_jj, local_matrix_in_out_row[jj]);
      }
    }
    for (size_t ii = 0; ii < rows_out; ++ii) {
      const auto& local_matrix_out_in_row = local_matrix_out_in[ii];
      const auto& local_matrix_out_out_row = local_matrix_out_out[ii];
      const size_t global_ii = global_row_indices_out[ii];
      for (size_t jj = 0; jj < cols_in; ++jj) {
        const size_t global_jj = global_col_indices_in[jj];
        global_matrix_out_in.add_to_entry(global_ii, global_jj, local_matrix_out_in_row[jj]);
      }
      for (size_t jj = 0; jj < cols_out; ++jj) {
        const size_t global_jj = global_col_indices_out[jj];
        global_matrix_out_out.add_to_entry(global_ii, global_jj, local_matrix_out_out_row[jj]);
      }
    }
  } // ... assemble(...)

  void assemble(const TestSpaceType& inner_test_space,
                const AnsatzSpaceType& inner_ansatz_space,
                const OuterTestSpaceType& outer_test_space,
                const OuterAnsatzSpaceType& outer_ansatz_space,
                const IntersectionType& intersection,
                MatrixType& global_matrix) const
  {
    assemble(inner_test_space,
             inner_ansatz_space,
             outer_test_space,
             outer_ansatz_space,
             intersection,
             global_matrix,
             global_matrix,
             global_matrix,
             global_matrix);
  }

private:
  const LocalCouplingTwoFormType& local_coupling_twoform_;
}; // class LocalCouplingTwoFormAssembler


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalCouplingOperatorApplicator : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_local_coupling_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalCouplingOperatorApplicator(const GridLayerType& grid_layer,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where)
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_layer, intersection);
  }

  virtual void
  apply_local(const IntersectionType& intersection, const EntityType& inside_entity, const EntityType& outside_entity)
  {
    local_operator_.apply(source_,
                          intersection,
                          *range_.local_discrete_function(inside_entity),
                          *range_.local_discrete_function(outside_entity));
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where_;
}; // class LocalCouplingOperatorApplicator


/**
 * \todo \attention Rename LocalBoundaryTwoFormAssemblerFunctor -> LocalBoundaryTwoFormAssembler after removing this
 *                  class!
 */
template <class TestSpace, class Intersection, class Matrix, class AnsatzSpace = TestSpace>
class DUNE_DEPRECATED_MSG("Use LocalBoundaryTwoFormAssemblerFunctor instead (13.05.2017)!")
    LocalBoundaryTwoFormAssembler
{
  static_assert(is_space<TestSpace>::value, "");
  static_assert(is_space<AnsatzSpace>::value, "");
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");

public:
  typedef TestSpace TestSpaceType;
  typedef AnsatzSpace AnsatzSpaceType;
  typedef Matrix MatrixType;
  typedef Intersection IntersectionType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalBoundaryTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        IntersectionType,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        FieldType>
      LocalBoundaryTwoFormType;

public:
  explicit LocalBoundaryTwoFormAssembler(const LocalBoundaryTwoFormType& local_twoform)
    : local_twoform_(local_twoform)
  {
  }

  void assemble(const TestSpaceType& test_space,
                const AnsatzSpaceType& ansatz_space,
                const IntersectionType& intersection,
                MatrixType& global_matrix) const
  {
    const auto entity = intersection.inside();
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    Dune::DynamicMatrix<FieldType> local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_twoform_.apply2(test_base, ansatz_base, intersection, local_matrix);
    // write local matrix to global
    // \todo: make mutable member, after SMP refactor
    const auto global_row_indices = test_space.mapper().globalIndices(entity);
    const auto global_col_indices = ansatz_space.mapper().globalIndices(entity);
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


template <class GridLayerType, class LocalOperatorType, class SourceType, class RangeType>
class LocalBoundaryOperatorApplicator : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_local_boundary_operator<LocalOperatorType>::value,
                "LocalOperatorType has to be derived from LocalCouplingOperatorInterface!");
  static_assert(XT::Functions::is_localizable_function<SourceType>::value,
                "SourceType has to be derived from XT::Functions::LocalizableFunctionInterface!");
  static_assert(is_discrete_function<RangeType>::value, "RangeType has to be a DiscreteFunctionDefault!");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;

  LocalBoundaryOperatorApplicator(const GridLayerType& grid_layer,
                                  const LocalOperatorType& local_operator,
                                  const SourceType& source,
                                  RangeType& range,
                                  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where)
    : grid_layer_(grid_layer)
    , local_operator_(local_operator)
    , source_(source)
    , range_(range)
    , where_(where)
  {
  }

  virtual bool apply_on(const GridLayerType& grid_layer, const IntersectionType& intersection) const
  {
    return where_.apply_on(grid_layer, intersection);
  }

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& /*outside_entity*/)
  {
    local_operator_.apply(source_, intersection, *range_.local_discrete_function(inside_entity));
  }

private:
  const GridLayerType& grid_layer_;
  const LocalOperatorType& local_operator_;
  const SourceType& source_;
  RangeType& range_;
  const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>& where_;
}; // class LocalBoundaryOperatorApplicator


/**
 * \todo \attention Rename LocalVolumeFunctionalAssemblerFunctor -> LocalVolumeFunctionalAssembler after removing this
 * class!
 */
template <class TestSpace, class Vector>
class DUNE_DEPRECATED_MSG("Use LocalFunctionalAssemblerFunctor instead (08.06.2017)!") LocalVolumeFunctionalAssembler
{
  static_assert(is_space<TestSpace>::value, "");
  static_assert(XT::LA::is_vector<Vector>::value, "");

public:
  typedef TestSpace TestSpaceType;
  typedef Vector VectorType;
  typedef typename TestSpaceType::EntityType EntityType;
  typedef typename VectorType::ScalarType FieldType;
  typedef LocalVolumeFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, FieldType>
      LocalVolumeFunctionalType;

  explicit LocalVolumeFunctionalAssembler(const LocalVolumeFunctionalType& local_volume_functional)
    : local_volume_functional_(local_volume_functional)
  {
  }

  void assemble(const TestSpaceType& test_space, const EntityType& entity, VectorType& global_vector) const
  {
    // prepare
    const size_t size = test_space.mapper().numDofs(entity);
    DynamicVector<FieldType> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_volume_functional_.apply(test_basis, local_vector);
    // write local vector to global
    // \todo: make mutable member, after SMP refactor
    const auto global_indices = test_space.mapper().globalIndices(entity);
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalVolumeFunctionalType& local_volume_functional_;
}; // class LocalVolumeFunctionalAssembler

/**
 * \todo \attention Rename LocalFaceFunctionalAssemblerFunctor -> LocalFaceFunctionalAssembler after removing this
 * class!
 */
template <class TestSpace, class Intersection, class Vector>
class DUNE_DEPRECATED_MSG("Use LocalFaceFunctionalAssemblerFunctor instead (26.06.2017)!") LocalFaceFunctionalAssembler
{
  static_assert(is_space<TestSpace>::value, "");
  static_assert(XT::Grid::is_intersection<Intersection>::value, "");
  static_assert(XT::LA::is_vector<Vector>::value, "");

public:
  typedef TestSpace TestSpaceType;
  typedef Intersection IntersectionType;
  typedef Vector VectorType;
  typedef typename VectorType::ScalarType FieldType;
  typedef LocalFaceFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, IntersectionType, FieldType>
      LocalFaceFunctionalType;

  explicit LocalFaceFunctionalAssembler(const LocalFaceFunctionalType& local_face_functional)
    : local_face_functional_(local_face_functional)
  {
  }

  void assemble(const TestSpaceType& test_space, const IntersectionType& intersection, VectorType& global_vector) const
  {
    // prepare
    const auto entity = intersection.inside();
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector<FieldType> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_face_functional_.apply(test_basis, intersection, local_vector);
    // write local vector to global
    // \todo: make mutable member, after SMP refactor
    const auto global_indices = test_space.mapper().globalIndices(entity);
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

private:
  const LocalFaceFunctionalType& local_face_functional_;
}; // class LocalFaceFunctionalAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMLBER_HH
