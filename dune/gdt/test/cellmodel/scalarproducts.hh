// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner (2019)

#ifndef DUNE_GDT_TEST_CELLMODEL_SCALARPRODUCTS_HH
#define DUNE_GDT_TEST_CELLMODEL_SCALARPRODUCTS_HH

#include <dune/istl/scalarproducts.hh>

using namespace Dune;

template <class VectorType, class MatrixType>
class MassMatrixScalarProduct : public ScalarProduct<VectorType>
{
public:
  //! export types
  typedef VectorType domain_type;
  typedef typename VectorType::field_type field_type;
  typedef typename FieldTraits<field_type>::real_type real_type;

  MassMatrixScalarProduct(const MatrixType& M)
    : M_(M)
    , tmp_vec_(M_.rows())
  {}

  /*! \brief Dot product of two vectors. In the complex case, the first argument is conjugated.
     It is assumed that the vectors are consistent on the interior+border
     partition.
   */
  virtual field_type dot(const VectorType& x, const VectorType& y) override final
  {
    M_.mv(x, tmp_vec_);
    return tmp_vec_.dot(y);
  }

  /*! \brief Norm of a right-hand side vector.
     The vector must be consistent on the interior+border partition
   */
  virtual real_type norm(const VectorType& x) override final
  {
    return std::sqrt(dot(x, x));
  }

  //! Category of the scalar product (see SolverCategory::Category)
  virtual SolverCategory::Category category() const override final
  {
    return SolverCategory::sequential;
  }

private:
  const MatrixType& M_;
  mutable VectorType tmp_vec_;
};


template <class VectorType, class MatrixType>
class PfieldScalarProduct : public ScalarProduct<VectorType>
{
public:
  //! export types
  typedef VectorType domain_type;
  typedef typename VectorType::field_type field_type;
  typedef typename FieldTraits<field_type>::real_type real_type;
  using ConstVectorViewType = XT::LA::ConstVectorView<VectorType>;

  PfieldScalarProduct(const MatrixType& M)
    : M_(M)
    , size_phi_(M_.rows())
    , tmp_vec_(size_phi_, 0., 0)
    , tmp_vec2_(size_phi_, 0., 0)
  {}

  virtual field_type dot(const VectorType& x, const VectorType& y) override final
  {
    const ConstVectorViewType y_phi(y, 0, size_phi_);
    const ConstVectorViewType y_phinat(y, size_phi_, 2 * size_phi_);
    const ConstVectorViewType y_mu(y, 2 * size_phi_, 3 * size_phi_);
    auto& x_phi = tmp_vec_;
    for (size_t ii = 0; ii < size_phi_; ++ii)
      x_phi[ii] = x[ii];
    M_.mv(x_phi, tmp_vec2_);
    field_type ret = y_phi.dot(tmp_vec2_);
    auto& x_phinat = tmp_vec_;
    for (size_t ii = 0; ii < size_phi_; ++ii)
      x_phinat[ii] = x[size_phi_ + ii];
    M_.mv(x_phinat, tmp_vec2_);
    ret += y_phinat.dot(tmp_vec2_);
    auto& x_mu = tmp_vec_;
    for (size_t ii = 0; ii < size_phi_; ++ii)
      x_mu[ii] = x[2 * size_phi_ + ii];
    M_.mv(x_mu, tmp_vec2_);
    ret += y_mu.dot(tmp_vec2_);
    return ret;
  }

  /*! \brief Norm of a right-hand side vector.
     The vector must be consistent on the interior+border partition
   */
  virtual real_type norm(const VectorType& x) override final
  {
    return std::sqrt(dot(x, x));
  }

  //! Category of the scalar product (see SolverCategory::Category)
  virtual SolverCategory::Category category() const override final
  {
    return SolverCategory::sequential;
  }

private:
  const MatrixType& M_;
  const size_t size_phi_;
  mutable VectorType tmp_vec_;
  mutable VectorType tmp_vec2_;
};

#endif // DUNE_GDT_TEST_CELLMODEL_SCALARPRODUCTS_HH
