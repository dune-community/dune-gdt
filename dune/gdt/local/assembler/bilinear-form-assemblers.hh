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

#include <dune/gdt/local/bilinear-forms/interfaces.hh>
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
class LocalElementBilinearFormAssembler : public XT::Grid::ElementFunctor<GridView>
{
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalElementBilinearFormAssembler<Matrix, GridView, t_r, t_rC, TR, TGV, AGV, a_r, a_rC, AR>;
  using BaseType = XT::Grid::ElementFunctor<GridView>;

public:
  using typename BaseType::ElementType;
  using MatrixType = Matrix;
  using FieldType = typename MatrixType::ScalarType;
  using TestSpaceType = SpaceInterface<TGV, t_r, t_rC, TR>;
  using AnsatzSpaceType = SpaceInterface<AGV, a_r, a_rC, AR>;
  using LocalBilinearFormType = LocalElementBilinearFormInterface<ElementType, t_r, t_rC, TR, FieldType, a_r, a_rC, AR>;

  LocalElementBilinearFormAssembler(const TestSpaceType& test_space,
                                    const AnsatzSpaceType& ansatz_space,
                                    const LocalBilinearFormType& local_two_form,
                                    MatrixType& global_matrix,
                                    const XT::Common::Parameter& param = {})
    : BaseType()
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , local_bilinear_form_(local_two_form.copy())
    , global_matrix_(global_matrix)
    , param_(param)
    , local_matrix_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , global_test_indices_(test_space_.mapper().max_local_size())
    , global_ansatz_indices_(ansatz_space_.mapper().max_local_size())
    , test_basis_(test_space_.basis().localize())
    , ansatz_basis_(ansatz_space_.basis().localize())
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

  LocalElementBilinearFormAssembler(const ThisType& other)
    : BaseType()
    , test_space_(other.test_space_)
    , ansatz_space_(other.ansatz_space_)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
    , global_matrix_(other.global_matrix_)
    , param_(other.param_)
    , local_matrix_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , global_test_indices_(test_space_.mapper().max_local_size())
    , global_ansatz_indices_(ansatz_space_.mapper().max_local_size())
    , test_basis_(test_space_.basis().localize())
    , ansatz_basis_(ansatz_space_.basis().localize())
  {
  }

  LocalElementBilinearFormAssembler(ThisType&& source) = default;

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const ElementType& element) override final
  {
    // apply bilinear form
    test_basis_->bind(element);
    ansatz_basis_->bind(element);
    local_bilinear_form_->apply2(*test_basis_, *ansatz_basis_, local_matrix_, param_);
    // copy local matrix to global matrix
    test_space_.mapper().global_indices(element, global_test_indices_);
    ansatz_space_.mapper().global_indices(element, global_ansatz_indices_);
    for (size_t ii = 0; ii < test_basis_->size(param_); ++ii)
      for (size_t jj = 0; jj < ansatz_basis_->size(param_); ++jj)
        global_matrix_.add_to_entry(global_test_indices_[ii], global_ansatz_indices_[jj], local_matrix_[ii][jj]);
  } // ... apply_local(...)

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  MatrixType& global_matrix_;
  XT::Common::Parameter param_;
  DynamicMatrix<FieldType> local_matrix_;
  DynamicVector<size_t> global_test_indices_;
  DynamicVector<size_t> global_ansatz_indices_;
  mutable std::unique_ptr<typename TestSpaceType::GlobalBasisType::LocalizedBasisType> test_basis_;
  mutable std::unique_ptr<typename AnsatzSpaceType::GlobalBasisType::LocalizedBasisType> ansatz_basis_;
}; // class LocalElementBilinearFormAssembler


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
class LocalIntersectionBilinearFormAssembler : public XT::Grid::IntersectionFunctor<GridView>
{
  static_assert(XT::LA::is_matrix<Matrix>::value, "");
  static_assert(XT::Grid::is_view<GridView>::value, "");

  using ThisType = LocalIntersectionBilinearFormAssembler<Matrix, GridView, t_r, t_rC, TR, TGV, AGV, a_r, a_rC, AR>;
  using BaseType = XT::Grid::IntersectionFunctor<GridView>;

public:
  using typename BaseType::ElementType;
  using typename BaseType::IntersectionType;
  using I = IntersectionType;
  using MatrixType = Matrix;
  using FieldType = typename MatrixType::ScalarType;
  using TestSpaceType = SpaceInterface<TGV, t_r, t_rC, TR>;
  using AnsatzSpaceType = SpaceInterface<AGV, a_r, a_rC, AR>;
  using LocalBilinearFormType = LocalIntersectionBilinearFormInterface<I, t_r, t_rC, TR, FieldType, a_r, a_rC, AR>;

  LocalIntersectionBilinearFormAssembler(const TestSpaceType& test_space,
                                         const AnsatzSpaceType& ansatz_space,
                                         const LocalBilinearFormType& local_two_form,
                                         MatrixType& global_matrix,
                                         const XT::Common::Parameter& param = {})
    : BaseType()
    , test_space_(test_space)
    , ansatz_space_(ansatz_space)
    , local_bilinear_form_(local_two_form.copy())
    , global_matrix_(global_matrix)
    , param_(param)
    , local_matrix_in_in_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_in_out_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_out_in_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_out_out_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , global_test_indices_in_(test_space_.mapper().max_local_size())
    , global_test_indices_out_(test_space_.mapper().max_local_size())
    , global_ansatz_indices_in_(ansatz_space_.mapper().max_local_size())
    , global_ansatz_indices_out_(ansatz_space_.mapper().max_local_size())
    , test_basis_inside_(test_space_.basis().localize())
    , test_basis_outside_(test_space_.basis().localize())
    , ansatz_basis_inside_(ansatz_space_.basis().localize())
    , ansatz_basis_outside_(ansatz_space_.basis().localize())
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

  LocalIntersectionBilinearFormAssembler(const ThisType& other)
    : BaseType()
    , test_space_(other.test_space_)
    , ansatz_space_(other.ansatz_space_)
    , local_bilinear_form_(other.local_bilinear_form_->copy())
    , global_matrix_(other.global_matrix_)
    , param_(other.param_)
    , local_matrix_in_in_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_in_out_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_out_in_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , local_matrix_out_out_(test_space_.mapper().max_local_size(), ansatz_space_.mapper().max_local_size())
    , global_test_indices_in_(test_space_.mapper().max_local_size())
    , global_test_indices_out_(test_space_.mapper().max_local_size())
    , global_ansatz_indices_in_(ansatz_space_.mapper().max_local_size())
    , global_ansatz_indices_out_(ansatz_space_.mapper().max_local_size())
    , test_basis_inside_(test_space_.basis().localize())
    , test_basis_outside_(test_space_.basis().localize())
    , ansatz_basis_inside_(ansatz_space_.basis().localize())
    , ansatz_basis_outside_(ansatz_space_.basis().localize())
  {
  }

  LocalIntersectionBilinearFormAssembler(ThisType&& source) = default;

  BaseType* copy() override final
  {
    return new ThisType(*this);
  }

  void apply_local(const IntersectionType& intersection,
                   const ElementType& inside_element,
                   const ElementType& outside_element) override final
  {
    // apply bilinear form
    test_basis_inside_->bind(inside_element);
    ansatz_basis_inside_->bind(inside_element);
    test_basis_outside_->bind(outside_element);
    ansatz_basis_outside_->bind(outside_element);
    local_bilinear_form_->apply2(intersection,
                                 *test_basis_inside_,
                                 *ansatz_basis_inside_,
                                 *test_basis_outside_,
                                 *ansatz_basis_outside_,
                                 local_matrix_in_in_,
                                 local_matrix_in_out_,
                                 local_matrix_out_in_,
                                 local_matrix_out_out_,
                                 param_);
    // copy local matrices to global matrix
    test_space_.mapper().global_indices(inside_element, global_test_indices_in_);
    test_space_.mapper().global_indices(outside_element, global_test_indices_out_);
    ansatz_space_.mapper().global_indices(inside_element, global_ansatz_indices_in_);
    ansatz_space_.mapper().global_indices(outside_element, global_ansatz_indices_out_);
    for (size_t ii = 0; ii < test_basis_inside_->size(param_); ++ii) {
      for (size_t jj = 0; jj < ansatz_basis_inside_->size(param_); ++jj)
        global_matrix_.add_to_entry(
            global_test_indices_in_[ii], global_ansatz_indices_in_[jj], local_matrix_in_in_[ii][jj]);
      for (size_t jj = 0; jj < ansatz_basis_outside_->size(param_); ++jj)
        global_matrix_.add_to_entry(
            global_test_indices_in_[ii], global_ansatz_indices_out_[jj], local_matrix_in_out_[ii][jj]);
    }
    for (size_t ii = 0; ii < test_basis_outside_->size(param_); ++ii) {
      for (size_t jj = 0; jj < ansatz_basis_inside_->size(param_); ++jj)
        global_matrix_.add_to_entry(
            global_test_indices_out_[ii], global_ansatz_indices_in_[jj], local_matrix_out_in_[ii][jj]);
      for (size_t jj = 0; jj < ansatz_basis_outside_->size(param_); ++jj)
        global_matrix_.add_to_entry(
            global_test_indices_out_[ii], global_ansatz_indices_out_[jj], local_matrix_out_out_[ii][jj]);
    }
  } // ... apply_local(...)

private:
  const TestSpaceType& test_space_;
  const AnsatzSpaceType& ansatz_space_;
  const std::unique_ptr<LocalBilinearFormType> local_bilinear_form_;
  MatrixType& global_matrix_;
  XT::Common::Parameter param_;
  DynamicMatrix<FieldType> local_matrix_in_in_;
  DynamicMatrix<FieldType> local_matrix_in_out_;
  DynamicMatrix<FieldType> local_matrix_out_in_;
  DynamicMatrix<FieldType> local_matrix_out_out_;
  DynamicVector<size_t> global_test_indices_in_;
  DynamicVector<size_t> global_test_indices_out_;
  DynamicVector<size_t> global_ansatz_indices_in_;
  DynamicVector<size_t> global_ansatz_indices_out_;
  mutable std::unique_ptr<typename TestSpaceType::GlobalBasisType::LocalizedBasisType> test_basis_inside_;
  mutable std::unique_ptr<typename TestSpaceType::GlobalBasisType::LocalizedBasisType> test_basis_outside_;
  mutable std::unique_ptr<typename AnsatzSpaceType::GlobalBasisType::LocalizedBasisType> ansatz_basis_inside_;
  mutable std::unique_ptr<typename AnsatzSpaceType::GlobalBasisType::LocalizedBasisType> ansatz_basis_outside_;
}; // class LocalIntersectionBilinearFormAssembler


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_ASSEMBLER_TWO_FORM_ASSEMBLERS_HH
