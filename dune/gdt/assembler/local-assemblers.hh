// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)
//   Rene Milk       (2018)
//   Tim Keil        (2017)

#ifndef DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH
#define DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/local/functionals/interfaces.hh>
#include <dune/gdt/type_traits.hh>

namespace Dune {
namespace GDT {


/**
 * \todo \attention Rename LocalVolumeTwoFormAssemblerFunctor -> LocalVolumeTwoFormAssembler after removing the latter!
 */
template <class TestSpaceType,
          class MatrixType,
          class GridLayerType = typename TestSpaceType::GridLayerType,
          class AnsatzSpaceType = TestSpaceType>
class LocalVolumeTwoFormAssemblerFunctor : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "");
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");
  static_assert(is_space<AnsatzSpaceType>::value, "");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalVolumeTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                      typename AnsatzSpaceType::BaseFunctionSetType,
                                      FieldType>
      LocalVolumeTwoFormType;

  static void assemble(const TestSpaceType& test_space,
                       const AnsatzSpaceType& ansatz_space,
                       const LocalVolumeTwoFormType& local_volume_two_form,
                       const EntityType& entity,
                       MatrixType& global_matrix)
  {

    DUNE_THROW_IF(global_matrix.rows() != test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix.rows() = " << global_matrix.rows() << "\n  "
                                            << "test_space.mapper().size()"
                                            << test_space.mapper().size());
    DUNE_THROW_IF(global_matrix.cols() != ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix.cols() = " << global_matrix.cols() << "\n  "
                                            << "ansatz_space.mapper().size()"
                                            << ansatz_space.mapper().size());
    // prepare
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    DynamicMatrix<FieldType> local_matrix(rows, cols, 0.); // \todo: make mutable member, after SMP refactor
    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_volume_two_form.apply2(test_base, ansatz_base, local_matrix);
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

  LocalVolumeTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                                     const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                     const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                                     const LocalVolumeTwoFormType& local_volume_two_form,
                                     MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_volume_two_form_(local_volume_two_form)
    , matrix_(matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    assemble(*test_space_, *ansatz_space_, local_volume_two_form_, entity, matrix_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  const LocalVolumeTwoFormType& local_volume_two_form_;
  MatrixType& matrix_;
}; // class LocalVolumeTwoFormAssembler


/**
 * \todo \attention Rename LocalCouplingTwoFormAssemblerFunctor -> LocalCouplingTwoFormAssembler after removing the
 *                  latter!
 */
template <class TestSpaceType,
          class MatrixType,
          class GridLayerType = typename TestSpaceType::GridLayerType,
          class AnsatzSpaceType = TestSpaceType,
          class OuterTestSpaceType = TestSpaceType,
          class OuterAnsatzSpaceType = AnsatzSpaceType>
class LocalCouplingTwoFormAssemblerFunctor : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "");
  static_assert(XT::LA::is_matrix<MatrixType>::value, "");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");
  static_assert(is_space<AnsatzSpaceType>::value, "");
  static_assert(is_space<OuterTestSpaceType>::value, "");
  static_assert(is_space<OuterAnsatzSpaceType>::value, "");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalCouplingTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        IntersectionType,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        typename OuterTestSpaceType::BaseFunctionSetType,
                                        typename OuterAnsatzSpaceType::BaseFunctionSetType,
                                        FieldType>
      LocalCouplingTwoFormType;

  static void assemble(const TestSpaceType& inner_test_space,
                       const AnsatzSpaceType& inner_ansatz_space,
                       const OuterTestSpaceType& outer_test_space,
                       const OuterAnsatzSpaceType& outer_ansatz_space,
                       const IntersectionType& intersection,
                       const LocalCouplingTwoFormType& local_coupling_two_form,
                       MatrixType& global_matrix_in_in,
                       MatrixType& global_matrix_out_out,
                       MatrixType& global_matrix_in_out,
                       MatrixType& global_matrix_out_in)
  {
    DUNE_THROW_IF(global_matrix_in_in.rows() != inner_test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_in_in.rows() = " << global_matrix_in_in.rows() << "\n  "
                                                  << "inner_test_space.mapper().size()"
                                                  << inner_test_space.mapper().size());
    DUNE_THROW_IF(global_matrix_in_in.cols() != inner_ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_in_in.cols() = " << global_matrix_in_in.cols() << "\n  "
                                                  << "inner_ansatz_space.mapper().size()"
                                                  << inner_ansatz_space.mapper().size());
    DUNE_THROW_IF(global_matrix_out_out.rows() != outer_test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_out_out.rows() = " << global_matrix_out_out.rows() << "\n  "
                                                    << "outer_test_space.mapper().size()"
                                                    << outer_test_space.mapper().size());
    DUNE_THROW_IF(global_matrix_out_out.cols() != outer_ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_out_out.cols() = " << global_matrix_out_out.cols() << "\n  "
                                                    << "outer_ansatz_space.mapper().size()"
                                                    << outer_ansatz_space.mapper().size());
    DUNE_THROW_IF(global_matrix_in_out.rows() != inner_test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_in_out.rows() = " << global_matrix_in_out.rows() << "\n  "
                                                   << "inner_test_space.mapper().size()"
                                                   << inner_test_space.mapper().size());
    DUNE_THROW_IF(global_matrix_in_out.cols() != outer_ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_in_out.cols() = " << global_matrix_in_out.cols() << "\n  "
                                                   << "outer_ansatz_space.mapper().size()"
                                                   << outer_ansatz_space.mapper().size());
    DUNE_THROW_IF(global_matrix_out_in.rows() != outer_test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_out_in.rows() = " << global_matrix_out_in.rows() << "\n  "
                                                   << "outer_test_space.mapper().size()"
                                                   << outer_test_space.mapper().size());
    DUNE_THROW_IF(global_matrix_out_in.cols() != inner_ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_out_in.cols() = " << global_matrix_out_in.cols() << "\n  "
                                                   << "inner_ansatz_space.mapper().size()"
                                                   << inner_ansatz_space.mapper().size());
    // prepare
    const auto entity = intersection.inside();
    const auto neighbor = intersection.outside();
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
    local_coupling_two_form.apply2(test_base_in,
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

  static void assemble(const TestSpaceType& inner_test_space,
                       const AnsatzSpaceType& inner_ansatz_space,
                       const OuterTestSpaceType& outer_test_space,
                       const OuterAnsatzSpaceType& outer_ansatz_space,
                       const IntersectionType& intersection,
                       const LocalCouplingTwoFormType& local_coupling_two_form,
                       MatrixType& global_matrix)
  {
    assemble(inner_test_space,
             inner_ansatz_space,
             outer_test_space,
             outer_ansatz_space,
             intersection,
             local_coupling_two_form,
             global_matrix,
             global_matrix,
             global_matrix,
             global_matrix);
  }

  template <typename TestSpace,
            typename AnsatzSpace,
            typename = typename std::enable_if<(std::is_same<TestSpace, OuterTestSpaceType>::value)
                                               && (std::is_same<AnsatzSpace, OuterAnsatzSpaceType>::value)
                                               && sizeof(TestSpace)
                                               && sizeof(AnsatzSpace)>::type>
  LocalCouplingTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpace>& test_space,
                                       const XT::Common::PerThreadValue<const AnsatzSpace>& ansatz_space,
                                       const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                       const LocalCouplingTwoFormType& local_coupling_two_form,
                                       MatrixType& in_in_matrix,
                                       MatrixType& out_out_matrix,
                                       MatrixType& in_out_matrix,
                                       MatrixType& out_in_matrix)
    : inner_test_space_(test_space)
    , inner_ansatz_space_(ansatz_space)
    , outer_test_space_(test_space)
    , outer_ansatz_space_(ansatz_space)
    , where_(where)
    , local_coupling_two_form_(local_coupling_two_form)
    , in_in_matrix_(in_in_matrix)
    , out_out_matrix_(out_out_matrix)
    , in_out_matrix_(in_out_matrix)
    , out_in_matrix_(out_in_matrix)
  {
  }

  template <typename TestSpace,
            typename AnsatzSpace,
            typename = typename std::enable_if<(std::is_same<TestSpace, OuterTestSpaceType>::value)
                                               && (std::is_same<AnsatzSpace, OuterAnsatzSpaceType>::value)
                                               && sizeof(TestSpace)
                                               && sizeof(AnsatzSpace)>::type>
  LocalCouplingTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpace>& test_space,
                                       const XT::Common::PerThreadValue<const AnsatzSpace>& ansatz_space,
                                       const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                       const LocalCouplingTwoFormType& local_coupling_two_form,
                                       MatrixType& matrix)
    : inner_test_space_(test_space)
    , inner_ansatz_space_(ansatz_space)
    , outer_test_space_(test_space)
    , outer_ansatz_space_(ansatz_space)
    , where_(where)
    , local_coupling_two_form_(local_coupling_two_form)
    , in_in_matrix_(matrix)
    , out_out_matrix_(matrix)
    , in_out_matrix_(matrix)
    , out_in_matrix_(matrix)
  {
  }

  LocalCouplingTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space,
                                       const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space,
                                       const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space,
                                       const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space,
                                       const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                       const LocalCouplingTwoFormType& local_coupling_two_form,
                                       MatrixType& in_in_matrix,
                                       MatrixType& out_out_matrix,
                                       MatrixType& in_out_matrix,
                                       MatrixType& out_in_matrix)
    : inner_test_space_(inner_test_space)
    , inner_ansatz_space_(inner_ansatz_space)
    , outer_test_space_(outer_test_space)
    , outer_ansatz_space_(outer_ansatz_space)
    , where_(where)
    , local_coupling_two_form_(local_coupling_two_form)
    , in_in_matrix_(in_in_matrix)
    , out_out_matrix_(out_out_matrix)
    , in_out_matrix_(in_out_matrix)
    , out_in_matrix_(out_in_matrix)
  {
  }

  LocalCouplingTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space,
                                       const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space,
                                       const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space,
                                       const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space,
                                       const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                       const LocalCouplingTwoFormType& local_coupling_two_form,
                                       MatrixType& matrix)
    : inner_test_space_(inner_test_space)
    , inner_ansatz_space_(inner_ansatz_space)
    , outer_test_space_(outer_test_space)
    , outer_ansatz_space_(outer_ansatz_space)
    , where_(where)
    , local_coupling_two_form_(local_coupling_two_form)
    , in_in_matrix_(matrix)
    , out_out_matrix_(matrix)
    , in_out_matrix_(matrix)
    , out_in_matrix_(matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    assemble(*inner_test_space_,
             *inner_ansatz_space_,
             *outer_test_space_,
             *outer_ansatz_space_,
             intersection,
             local_coupling_two_form_,
             in_in_matrix_,
             out_out_matrix_,
             in_out_matrix_,
             out_in_matrix_);
  } // ... apply_local(...)

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& inner_test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& inner_ansatz_space_;
  const XT::Common::PerThreadValue<const OuterTestSpaceType>& outer_test_space_;
  const XT::Common::PerThreadValue<const OuterAnsatzSpaceType>& outer_ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalCouplingTwoFormType& local_coupling_two_form_;
  MatrixType& in_in_matrix_;
  MatrixType& out_out_matrix_;
  MatrixType& in_out_matrix_;
  MatrixType& out_in_matrix_;
}; // class LocalCouplingTwoFormAssemblerFunctor


/**
 * \todo \attention Rename LocalBoundaryTwoFormAssemblerFunctor -> LocalBoundaryTwoFormAssembler after removing the
 *                  latter!
 */
template <class TestSpaceType,
          class MatrixType,
          class GridLayerType = typename TestSpaceType::GridLayerType,
          class AnsatzSpaceType = TestSpaceType>
class LocalBoundaryTwoFormAssemblerFunctor : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef typename MatrixType::ScalarType FieldType;
  typedef LocalBoundaryTwoFormInterface<typename TestSpaceType::BaseFunctionSetType,
                                        IntersectionType,
                                        typename AnsatzSpaceType::BaseFunctionSetType,
                                        FieldType>
      LocalBoundaryTwoFormType;

  static void assemble(const TestSpaceType& test_space,
                       const AnsatzSpaceType& ansatz_space,
                       const IntersectionType& intersection,
                       const LocalBoundaryTwoFormType& local_boundary_two_form,
                       MatrixType& global_matrix)
  {
    DUNE_THROW_IF(global_matrix.rows() != test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix.rows() = " << global_matrix.rows() << "\n  "
                                            << "test_space.mapper().size()"
                                            << test_space.mapper().size());
    DUNE_THROW_IF(global_matrix.cols() != ansatz_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix.cols() = " << global_matrix.cols() << "\n  "
                                            << "ansatz_space.mapper().size()"
                                            << ansatz_space.mapper().size());
    // prepare
    const auto entity = intersection.inside();
    const size_t rows = test_space.mapper().numDofs(entity);
    const size_t cols = ansatz_space.mapper().numDofs(entity);
    // \todo: make mutable member, after SMP refactor
    Dune::DynamicMatrix<FieldType> local_matrix(rows, cols, 0.);
    // apply local two-form
    const auto test_base = test_space.base_function_set(entity);
    const auto ansatz_base = ansatz_space.base_function_set(entity);
    assert(test_base.size() == rows);
    assert(ansatz_base.size() == cols);
    local_boundary_two_form.apply2(test_base, ansatz_base, intersection, local_matrix);
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

  LocalBoundaryTwoFormAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& test_space,
                                       const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space,
                                       const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                       const LocalBoundaryTwoFormType& local_boundary_two_form,
                                       MatrixType& matrix)
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , where_(where)
    , local_boundary_two_form_(local_boundary_two_form)
    , matrix_(matrix)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    assemble(*test_space_, *ansatz_space_, intersection, local_boundary_two_form_, matrix_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& test_space_;
  const XT::Common::PerThreadValue<const AnsatzSpaceType>& ansatz_space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalBoundaryTwoFormType& local_boundary_two_form_;
  MatrixType& matrix_;
}; // class LocalBoundaryTwoFormAssemblerFunctor


/**
 * \todo \attention Rename LocalVolumenFunctionalAssemblerFunctor -> LocalFunctionalAssembler after removing the
 *                  latter!
 */
template <class TestSpaceType, class VectorType, class GridLayerType = typename TestSpaceType::GridLayerType>
class LocalVolumeFunctionalAssemblerFunctor : public XT::Grid::internal::Codim0Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "");
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");
  typedef XT::Grid::internal::Codim0Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  typedef typename VectorType::ScalarType FieldType;
  typedef LocalVolumeFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, FieldType>
      LocalVolumeFunctionalType;

  static void assemble(const TestSpaceType& test_space,
                       const LocalVolumeFunctionalType& local_volume_functional,
                       const EntityType& entity,
                       VectorType& global_vector)
  {
    DUNE_THROW_IF(global_vector.size() != test_space.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_vector.size() = " << global_vector.size() << "\n  "
                                            << "test_space.mapper().size()"
                                            << test_space.mapper().size());
    // prepare
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector<FieldType> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_volume_functional.apply(test_basis, local_vector);
    // write local vector to global
    const auto global_indices =
        test_space.mapper().globalIndices(entity); // \todo: make mutable member, after SMP refactor
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

  LocalVolumeFunctionalAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& space,
                                        const XT::Grid::ApplyOn::WhichEntity<GridLayerType>* where,
                                        const LocalVolumeFunctionalType& local_volume_functional,
                                        VectorType& vector)
    : space_(space)
    , where_(where)
    , local_volume_functional_(local_volume_functional)
    , vector_(vector)
  {
  }

  bool apply_on(const GridLayerType& gv, const EntityType& entity) const override final
  {
    return where_->apply_on(gv, entity);
  }

  void apply_local(const EntityType& entity) override final
  {
    assemble(*space_, local_volume_functional_, entity, vector_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichEntity<GridLayerType>> where_;
  const LocalVolumeFunctionalType& local_volume_functional_;
  VectorType& vector_;
}; // class LocalVolumeFunctionalAssemblerFunctor


/**
 * \todo \attention Rename LocalFaceFunctionalAssemblerFunctor -> LocalFaceFunctionalAssembler after removing the
 *                  latter!
 */
template <class TestSpaceType, class VectorType, class GridLayerType = typename TestSpaceType::GridLayerType>
class LocalFaceFunctionalAssemblerFunctor : public XT::Grid::internal::Codim1Object<GridLayerType>
{
  static_assert(is_space<TestSpaceType>::value, "");
  static_assert(XT::LA::is_vector<VectorType>::value, "");
  static_assert(XT::Grid::is_layer<GridLayerType>::value, "");
  typedef XT::Grid::internal::Codim1Object<GridLayerType> BaseType;

public:
  using typename BaseType::EntityType;
  using typename BaseType::IntersectionType;
  typedef typename VectorType::ScalarType FieldType;
  typedef LocalFaceFunctionalInterface<typename TestSpaceType::BaseFunctionSetType, IntersectionType, FieldType>
      LocalFaceFunctionalType;

  static void assemble(const TestSpaceType& test_space,
                       const IntersectionType& intersection,
                       const LocalFaceFunctionalType& local_face_functional,
                       VectorType& global_vector)
  {
    // \todo: if statement for global indices?
    // prepare
    const auto entity = intersection.inside();
    const size_t size = test_space.mapper().numDofs(entity);
    Dune::DynamicVector<FieldType> local_vector(size, 0.); // \todo: make mutable member, after SMP refactor
    // apply local functional
    const auto test_basis = test_space.base_function_set(entity);
    assert(test_basis.size() == size);
    local_face_functional.apply(test_basis, intersection, local_vector);
    // write local vector to global
    // \todo: make mutable member, after SMP refactor
    const auto global_indices = test_space.mapper().globalIndices(entity);
    assert(global_indices.size() == size);
    for (size_t jj = 0; jj < size; ++jj)
      global_vector.add_to_entry(global_indices[jj], local_vector[jj]);
  } // ... assemble(...)

  LocalFaceFunctionalAssemblerFunctor(const XT::Common::PerThreadValue<const TestSpaceType>& space,
                                      const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>* where,
                                      const LocalFaceFunctionalType& local_face_functional,
                                      VectorType& vector)
    : space_(space)
    , where_(where)
    , local_face_functional_(local_face_functional)
    , vector_(vector)
  {
  }

  bool apply_on(const GridLayerType& gv, const IntersectionType& intersection) const override final
  {
    return where_->apply_on(gv, intersection);
  }

  void apply_local(const IntersectionType& intersection,
                   const EntityType& /*inside_entity*/,
                   const EntityType& /*outside_entity*/) override final
  {
    assemble(*space_, intersection, local_face_functional_, vector_);
  }

private:
  const XT::Common::PerThreadValue<const TestSpaceType>& space_;
  const std::unique_ptr<const XT::Grid::ApplyOn::WhichIntersection<GridLayerType>> where_;
  const LocalFaceFunctionalType& local_face_functional_;
  VectorType& vector_;
}; // class LocalFaceFunctionalAssemblerFunctor


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH
