// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2017)

#ifndef DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH
#define DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH

#include <dune/xt/common/exceptions.hh>
#include <dune/xt/la/type_traits.hh>
#include <dune/xt/grid/walker/wrapper.hh>
#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/operators/interfaces.hh>
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
    if (global_matrix.rows() != test_space.mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
                 "global_matrix.rows() = " << global_matrix.rows() << "\n  "
                                           << "test_space.mapper().size()"
                                           << test_space.mapper().size());
    if (global_matrix.cols() != ansatz_space.mapper().size())
      DUNE_THROW(XT::Common::Exceptions::shapes_do_not_match,
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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_ASSEMBLER_LOCAL_ASSEMBLERS_HH
