// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_HH

#include <dune/stuff/common/disable_warnings.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>
#include <dune/stuff/common/reenable_warnings.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/crtp.hh>

namespace Dune {
namespace GDT {
namespace Spaces {


template <class ValueImp>
class ConstraintsInterface
{
public:
  typedef ValueImp ValueType;

  virtual ~ConstraintsInterface()
  {
  }

  virtual size_t rows() const = 0;

  virtual size_t cols() const = 0;

  virtual size_t global_row(const size_t ii) const = 0;

  virtual size_t global_col(const size_t jj) const = 0;

  virtual ValueType value(const size_t ii, const size_t jj) const = 0;
}; // class ConstraintsInterface


namespace Constraints {
namespace internal {


template <class ValueImp>
class Default : public ConstraintsInterface<ValueImp>
{
  typedef ConstraintsInterface<ValueImp> BaseType;

public:
  typedef typename BaseType::ValueType ValueType;

private:
  typedef DynamicVector<size_t> IndicesType;
  typedef DynamicMatrix<ValueType> ValuesType;

public:
  Default(const size_t rws, const size_t cls)
    : rows_(rws)
    , cols_(cls)
    , global_rows_(rows_)
    , global_cols_(cols_)
    , values_(rows_, cols_)
  {
  }

  virtual ~Default()
  {
  }

  virtual size_t rows() const DS_OVERRIDE
  {
    return rows_;
  }

  virtual size_t cols() const DS_OVERRIDE
  {
    return cols_;
  }

  void set_size(const size_t rr, const size_t cc)
  {
    rows_        = rr;
    cols_        = cc;
    bool changed = false;
    if (rows_ > global_rows_.size()) {
      global_rows_.resize(rows_, 0);
      changed = true;
    }
    if (cols_ > global_cols_.size()) {
      global_cols_.resize(cols_, 0);
      changed = true;
    }
    if (changed)
      values_.resize(global_rows_.size(), global_cols_.size(), 0);
    assert(global_rows_.size() == rows_);
    assert(global_cols_.size() == cols_);
  } // ... set_size(...)

  size_t& global_row(const size_t ii)
  {
    assert(ii < rows_);
    return global_rows_[ii];
  }

  virtual size_t global_row(const size_t ii) const DS_OVERRIDE
  {
    assert(ii < rows_);
    return global_rows_[ii];
  }

  size_t& global_col(const size_t jj)
  {
    assert(jj < cols_);
    return global_cols_[jj];
  }

  virtual size_t global_col(const size_t jj) const
  {
    assert(jj < cols_);
    return global_cols_[jj];
  }

  ValueType& value(const size_t ii, const size_t jj)
  {
    assert(ii < rows_);
    assert(jj < cols_);
    return values_[ii][jj];
  }

  virtual ValueType value(const size_t ii, const size_t jj) const DS_OVERRIDE
  {
    assert(ii < rows_);
    assert(jj < cols_);
    return values_[ii][jj];
  }

private:
  size_t rows_;
  size_t cols_;
  IndicesType global_rows_;
  IndicesType global_cols_;
  ValuesType values_;
}; // class Default


} // namespace internal


template <class IntersectionType, class ValueImp = double, bool setRow = true>
class Dirichlet : public internal::Default<ValueImp>
{
public:
  typedef internal::Default<ValueImp> BaseType;
  typedef Stuff::Grid::BoundaryInfoInterface<IntersectionType> BoundaryInfoType;

  Dirichlet(const BoundaryInfoType& bound_info, const size_t rws, const size_t cls)
    : BaseType(rws, cls)
    , boundary_info_(bound_info)
  {
  }

  virtual ~Dirichlet()
  {
  }

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

private:
  const BoundaryInfoType& boundary_info_;
}; // class Dirichlet


} // namespace Constraints
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
