// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_ASSEMBLER_TWO_FORM_ASSEMBLERS_HH
#define DUNE_GDT_LOCAL_ASSEMBLER_TWO_FORM_ASSEMBLERS_HH

#include <dune/xt/grid/functors/interfaces.hh>
#include <dune/xt/la/container/matrix-interface.hh>

#include <dune/gdt/local/operators/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


template <class Matrix,
          class GridView,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class TGV = GridView,
          class AGV = GridView,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class LocalElementTwoFormAssembler : public XT::Grid::ElementFunctor<GridView>
{
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalElementTwoFormAssembler<Matrix, GridView, t_r, t_rC, TR, TGV, AGV, a_r, a_rC, AR>;
  using BaseType = XT::Grid::ElementFunctor<GridView>;

public:
  using typename BaseType::ElementType;
  using MatrixType = Matrix;
  using FieldType = typename MatrixType::ScalarType;
  using TestSpaceType = SpaceInterface<TGV, t_r, t_rC, TR>;
  using AnsatzSpaceType = SpaceInterface<AGV, a_r, a_rC, AR>;
  using LocalTwoFormType = LocalElementTwoFormInterface<ElementType, t_r, t_rC, TR, FieldType, a_r, a_rC, AR>;

  LocalElementTwoFormAssembler(const TestSpaceType& test_space,
                               const AnsatzSpaceType& ansatz_space,
                               const LocalTwoFormType& local_two_form,
                               MatrixType& global_matrix,
                               const XT::Common::Parameter& param = {})
    : test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , local_two_form_(local_two_form.copy())
    , global_matrix_(global_matrix)
    , param_(param)
    , local_matrix_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , global_test_indices_(test_space_.mapper().max_local_size())
    , global_ansatz_indices_(ansatz_space_.mapper().max_local_size())
  {
    DUNE_THROW_IF(global_matrix_.rows() != test_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_.rows() = " << global_matrix_.rows() << "\n  "
                                             << "test_space_.mapper().size()"
                                             << test_space_.mapper().size());
    DUNE_THROW_IF(global_matrix_.cols() != ansatz_space_.mapper().size(),
                  XT::Common::Exceptions::shapes_do_not_match,
                  "global_matrix_.cols() = " << global_matrix_.cols() << "\n  "
                                             << "ansatz_space_.mapper().size()"
                                             << ansatz_space_.mapper().size());
  }

  LocalElementTwoFormAssembler(const ThisType& other) = default;
  LocalElementTwoFormAssembler(ThisType&& source) = default;

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const ElementType& element) override final
  {
    // apply two-form
    const auto test_basis = test_space_.basis().localize(element);
    const auto ansatz_basis = ansatz_space_.basis().localize(element);
    local_two_form_->apply2(*test_basis, *ansatz_basis, local_matrix_, param_);
    // copy local vector to global
    test_space_.mapper().global_indices(element, global_test_indices_);
    ansatz_space_.mapper().global_indices(element, global_ansatz_indices_);
    for (size_t ii = 0; ii < test_basis->size(param_); ++ii)
      for (size_t jj = 0; jj < ansatz_basis->size(param_); ++jj)
        global_matrix_.add_to_entry(global_test_indices_[ii], global_ansatz_indices_[jj], local_matrix_[ii][jj]);
  } // ... apply_local(...)

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
  const std::unique_ptr<LocalTwoFormType> local_two_form_;
  MatrixType& global_matrix_;
  XT::Common::Parameter param_;
  DynamicMatrix<FieldType> local_matrix_;
  DynamicVector<size_t> global_test_indices_;
  DynamicVector<size_t> global_ansatz_indices_;
}; // class LocalElementTwoFormAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_TWO_FORM_ASSEMBLERS_HH
