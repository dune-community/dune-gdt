// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_BILINEAR_FORMS_GENERIC_HH
#define DUNE_GDT_LOCAL_BILINEAR_FORMS_GENERIC_HH

#include <functional>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * See also LocalElementBilinearFormInterface for a description of the template arguments.
 *
 * \sa LocalElementBilinearFormInterface
 */
template <class E,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class GenericLocalElementBilinearForm : public LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = GenericLocalElementBilinearForm<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;
  using BaseType = LocalElementBilinearFormInterface<E, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericFunctionType = std::function<void(const LocalTestBasisType& /*test_basis*/,
                                                 const LocalAnsatzBasisType& /*ansatz_basis*/,
                                                 DynamicMatrix<F>& /*result*/,
                                                 const XT::Common::Parameter& /*param*/)>;

  GenericLocalElementBilinearForm(GenericFunctionType func, const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , func_(func)
  {}

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const LocalTestBasisType& test_basis,
              const LocalAnsatzBasisType& ansatz_basis,
              DynamicMatrix<F>& result,
              const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t rows = test_basis.size(param);
    const size_t cols = ansatz_basis.size(param);
    if (result.rows() < rows || result.cols() < cols)
      result.resize(rows, cols);
    result *= 0;
    // compute
    func_(test_basis, ansatz_basis, result, this->parse_parameter(param));
  } // ... apply2(...)

private:
  const GenericFunctionType func_;
}; // class GenericLocalElementBilinearForm


/**
 * See also LocalIntersectionBilinearFormInterface for a description of the template arguments.
 *
 * \sa LocalIntersectionBilinearFormInterface
 */
template <class I,
          size_t t_r = 1,
          size_t t_rC = 1,
          class TR = double,
          class F = double,
          size_t a_r = t_r,
          size_t a_rC = t_rC,
          class AR = TR>
class GenericLocalIntersectionBilinearForm
  : public LocalIntersectionBilinearFormInterface<I, t_r, t_rC, TR, F, a_r, a_rC, AR>
{
  using ThisType = GenericLocalIntersectionBilinearForm<I, t_r, t_rC, TR, F, a_r, a_rC, AR>;
  using BaseType = LocalIntersectionBilinearFormInterface<I, t_r, t_rC, TR, F, a_r, a_rC, AR>;

public:
  using typename BaseType::IntersectionType;
  using typename BaseType::LocalAnsatzBasisType;
  using typename BaseType::LocalTestBasisType;

  using GenericFunctionType = std::function<void(const IntersectionType& /*intersection*/,
                                                 const LocalTestBasisType& /*test_basis_inside*/,
                                                 const LocalAnsatzBasisType& /*ansatz_basis_inside*/,
                                                 const LocalTestBasisType& /*test_basis_outside*/,
                                                 const LocalAnsatzBasisType& /*ansatz_basis_outside*/,
                                                 DynamicMatrix<F>& /*result_in_in*/,
                                                 DynamicMatrix<F>& /*result_in_out*/,
                                                 DynamicMatrix<F>& /*result_out_in*/,
                                                 DynamicMatrix<F>& /*result_out_out*/,
                                                 const XT::Common::Parameter& /*param*/)>;

  GenericLocalIntersectionBilinearForm(GenericFunctionType func, const XT::Common::ParameterType& param_type = {})
    : BaseType(param_type)
    , func_(func)
  {}

  GenericLocalIntersectionBilinearForm(ThisType&& source) = default;

  std::unique_ptr<BaseType> copy() const override final
  {
    return std::make_unique<ThisType>(*this);
  }

  using BaseType::apply2;

  void apply2(const IntersectionType& intersection,
              const LocalTestBasisType& test_basis_inside,
              const LocalAnsatzBasisType& ansatz_basis_inside,
              const LocalTestBasisType& test_basis_outside,
              const LocalAnsatzBasisType& ansatz_basis_outside,
              DynamicMatrix<F>& result_in_in,
              DynamicMatrix<F>& result_in_out,
              DynamicMatrix<F>& result_out_in,
              DynamicMatrix<F>& result_out_out,
              const XT::Common::Parameter& param = {}) const override final
  {
    // prepare storage
    const size_t rows_in = test_basis_inside.size(param);
    const size_t rows_out = test_basis_outside.size(param);
    const size_t cols_in = ansatz_basis_inside.size(param);
    const size_t cols_out = ansatz_basis_outside.size(param);
    const auto ensure_size_and_clear = [](auto& m, const auto& r, const auto& c) {
      if (m.rows() < r || m.cols() < c)
        m.resize(r, c);
      m *= 0;
    };
    ensure_size_and_clear(result_in_in, rows_in, cols_in);
    ensure_size_and_clear(result_in_out, rows_in, cols_out);
    ensure_size_and_clear(result_out_in, rows_out, cols_in);
    ensure_size_and_clear(result_out_out, rows_out, cols_out);
    // compute
    func_(intersection,
          test_basis_inside,
          ansatz_basis_inside,
          test_basis_outside,
          ansatz_basis_outside,
          result_in_in,
          result_in_out,
          result_out_in,
          result_out_out,
          param);
  } // ... apply2(...)

private:
  const GenericFunctionType func_;
}; // class GenericLocalIntersectionBilinearForm


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_BILINEAR_FORMS_GENERIC_HH
