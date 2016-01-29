// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_SPACES_CONSTRAINTS_HH
#define DUNE_GDT_SPACES_CONSTRAINTS_HH

#include <ostream>

#include <dune/common/dynvector.hh>
#include <dune/common/dynmatrix.hh>

#include <dune/stuff/grid/boundaryinfo.hh>
#include <dune/stuff/common/crtp.hh>

namespace Dune {
namespace GDT {
namespace Spaces {


template< class Traits, class ValueImp = double >
class ConstraintsInterface
  : public Stuff::CRTPInterface< ConstraintsInterface< Traits, ValueImp >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef ValueImp                      ValueType;

  inline size_t rows() const
  {
    CHECK_CRTP(this->as_imp().rows());
    return this->as_imp().rows();
  }

  inline size_t cols() const
  {
    CHECK_CRTP(this->as_imp().cols());
    return this->as_imp().cols();
  }

  inline size_t global_row(const size_t ii) const
  {
    CHECK_CRTP(this->as_imp().global_row(ii));
    return this->as_imp().global_row(ii);
  }

  inline size_t global_col(const size_t jj) const
  {
    CHECK_CRTP(this->as_imp().global_col(jj));
    return this->as_imp().global_col(jj);
  }

  inline ValueType value(const size_t ii, const size_t jj) const
  {
    CHECK_CRTP(this->as_imp().value(ii, jj));
    return this->as_imp().value(ii, jj);
  }

private:
  template< class T, class V >
  friend std::ostream& operator<<(std::ostream& /*out*/, const ConstraintsInterface< T, V >& /*constraints*/);
}; // class ConstraintsInterface


template< class T, class V >
std::ostream& operator<<(std::ostream& out, const ConstraintsInterface< T, V >& constraints)
{
  out << "gdt.spaces.constraints of size " << constraints.rows() << " x " << constraints.cols() << ", containing\n";
  for (size_t ii = 0; ii < constraints.rows(); ++ii)
    for (size_t jj = 0; jj < constraints.cols(); ++jj)
      out << "  " << constraints.global_row(ii) << ", " << constraints.global_col(jj) << ": "
          << constraints.value(ii, jj) << "\n";
  return out;
} // ... operator<<(...)


namespace Constraints {
namespace internal {


template< class DerivedTraits, class ValueImp >
class Default
  : public ConstraintsInterface< DerivedTraits, ValueImp >
{
  typedef ConstraintsInterface< DerivedTraits, ValueImp > BaseType;
public:
  typedef DerivedTraits                Traits;
  typedef typename BaseType::ValueType ValueType;

private:
  typedef DynamicVector< size_t >    IndicesType;
  typedef DynamicMatrix< ValueType > ValuesType;

public:
  Default(const size_t rws, const size_t cls)
    : rows_(rws)
    , cols_(cls)
    , global_rows_(rows_)
    , global_cols_(cols_)
    , values_(rows_, cols_)
  {}

  size_t rows() const
  {
    return rows_;
  }

  size_t cols() const
  {
    return cols_;
  }

  void set_size(const size_t rr, const size_t cc)
  {
    rows_ = rr;
    cols_ = cc;
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
  } // ... set_size(...)

  size_t& global_row(const size_t ii)
  {
    assert(ii < std::min(rows_, global_rows_.size()));
    return global_rows_[ii];
  }

  size_t global_row(const size_t ii) const
  {
    assert(ii < std::min(rows_, global_rows_.size()));
    return global_rows_[ii];
  }

  size_t& global_col(const size_t jj)
  {
    assert(jj < std::min(cols_, global_cols_.size()));
    return global_cols_[jj];
  }

  size_t global_col(const size_t jj) const
  {
    assert(jj < std::min(cols_, global_cols_.size()));
    return global_cols_[jj];
  }

  ValueType& value(const size_t ii, const size_t jj)
  {
    assert(ii < std::min(rows_, global_rows_.size()));
    assert(jj < std::min(cols_, global_cols_.size()));
    return values_[ii][jj];
  }

  ValueType value(const size_t ii, const size_t jj) const
  {
    assert(ii < std::min(rows_, global_rows_.size()));
    assert(jj < std::min(cols_, global_cols_.size()));
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


template< class IntersectionType, class ValueImp = double >
class Dirichlet;


namespace internal {


template< class IntersectionType, class ValueImp >
class DirichletTraits
{
public:
  typedef Dirichlet< IntersectionType, ValueImp > derived_type;
};


} // namespace internal


template< class IntersectionType, class ValueImp >
class Dirichlet
  : public internal::Default< internal::DirichletTraits< IntersectionType, ValueImp >, ValueImp >
{
  typedef internal::Default< internal::DirichletTraits< IntersectionType, ValueImp >, ValueImp > BaseType;
public:
  typedef Stuff::Grid::BoundaryInfoInterface< IntersectionType > BoundaryInfoType;

  Dirichlet(const BoundaryInfoType& bnd_info, const size_t rws, const size_t cls, const bool set_rw = true)
    : BaseType(rws, cls)
    , boundary_info_(bnd_info)
    , set_row_(set_rw)
  {}

  const BoundaryInfoType& boundary_info() const
  {
    return boundary_info_;
  }

  bool set_row() const
  {
    return set_row_;
  }

private:
  const BoundaryInfoType& boundary_info_;
  const bool set_row_;
}; // class Dirichlet


} // namespace Constraints
} // namespace Spaces
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_CONSTRAINTS_HH
