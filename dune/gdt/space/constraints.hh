// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACE_CONSTRAINTS_HH
#define DUNE_GDT_SPACE_CONSTRAINTS_HH

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/stuff/grid/boundaryinfo.hh>

namespace Dune {
namespace GDT {


namespace Constraints {


template <class RangeFieldImp = double>
class LocalDefault
{
public:
  typedef RangeFieldImp RangeFieldType;

private:
  typedef Dune::DynamicVector<size_t> IndicesType;
  typedef Dune::DynamicMatrix<RangeFieldType> ValuesType;

public:
  LocalDefault(const size_t numRows, const size_t numCols)
    : rows_(numRows)
    , cols_(numCols)
    , globalRows_(rows_)
    , globalCols_(cols_)
    , values_(rows_, cols_)
  {
  }

  virtual ~LocalDefault()
  {
  }

  virtual size_t rows() const
  {
    return rows_;
  }

  virtual size_t cols() const
  {
    return cols_;
  }

  void setSize(const size_t rr, const size_t cc)
  {
    rows_        = rr;
    cols_        = cc;
    bool changed = false;
    if (rows_ > globalRows_.size()) {
      globalRows_.resize(rows_, 0);
      changed = true;
    }
    if (cols_ > globalCols_.size()) {
      globalCols_.resize(cols_, 0);
      changed = true;
    }
    if (changed)
      values_.resize(globalRows_.size(), globalCols_.size(), 0);
  } // ... setSize(...)

  virtual size_t& globalRow(const size_t ii)
  {
    assert(ii < std::min(rows_, globalRows_.size()));
    return globalRows_[ii];
  }

  virtual const size_t& globalRow(const size_t ii) const
  {
    assert(ii < std::min(rows_, globalRows_.size()));
    return globalRows_[ii];
  }

  virtual size_t& globalCol(const size_t jj)
  {
    assert(jj < std::min(cols_, globalCols_.size()));
    return globalCols_[jj];
  }

  virtual const size_t& globalCol(const size_t jj) const
  {
    assert(jj < std::min(cols_, globalCols_.size()));
    return globalCols_[jj];
  }

  virtual RangeFieldType& value(const size_t ii, const size_t jj)
  {
    assert(ii < std::min(rows_, globalRows_.size()));
    assert(jj < std::min(cols_, globalCols_.size()));
    return values_[ii][jj];
  }

  virtual const RangeFieldType& value(const size_t ii, const size_t jj) const
  {
    assert(ii < std::min(rows_, globalRows_.size()));
    assert(jj < std::min(cols_, globalCols_.size()));
    return values_[ii][jj];
  }

private:
  size_t rows_;
  size_t cols_;
  IndicesType globalRows_;
  IndicesType globalCols_;
  ValuesType values_;
}; // class LocalDefault


template <class GridViewType, class RangeFieldImp = double, bool setRow = true>
class Dirichlet : public LocalDefault<RangeFieldImp>
{
public:
  typedef LocalDefault<RangeFieldImp> BaseType;
  typedef Dune::Stuff::GridboundaryInterface<GridViewType> GridBoundaryType;

  Dirichlet(const GridBoundaryType& gB, const size_t numRows, const size_t numCols)
    : BaseType(numRows, numCols)
    , gridBoundary_(gB)
  {
  }

  const GridBoundaryType& gridBoundary() const
  {
    return gridBoundary_;
  }

private:
  const GridBoundaryType& gridBoundary_;
}; // class Dirichlet


} // namespace Constraints
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACE_CONSTRAINTS_HH
