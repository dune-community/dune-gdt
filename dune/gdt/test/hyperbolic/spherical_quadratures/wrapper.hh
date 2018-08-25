// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Tobias Leibner  (2017)

#ifndef DUNE_GDT_TEST_HYPERBOLIC_SPHERICALQUADRATURES_WRAPPER_HH
#define DUNE_GDT_TEST_HYPERBOLIC_SPHERICALQUADRATURES_WRAPPER_HH

#include <numeric>
#include <vector>

#include <dune/geometry/quadraturerules.hh>

namespace Dune {
namespace GDT {
namespace Hyperbolic {
namespace Problems {


// forward
template <class FieldType, size_t dimDomain>
class QuadraturesWrapper;

// Wrapper to be able to iterate over all quadratures in QuadraturesWrapper with one loop
template <class FieldType, size_t dimDomain>
class MergedQuadrature
{
public:
  using QuadraturesWrapperType = QuadraturesWrapper<FieldType, dimDomain>;
  using QuadratureType = typename QuadraturesWrapperType::QuadratureType;
  using QuadPointType = typename QuadraturesWrapper<FieldType, dimDomain>::QuadPointType;

  MergedQuadrature(const QuadraturesWrapperType& quadratures)
    : quadratures_(quadratures)
  {
  }

  class MergedQuadratureIterator
  {
  public:
    MergedQuadratureIterator()
    {
    }

    MergedQuadratureIterator(const QuadraturesWrapperType& quadratures,
                             const size_t first_index,
                             const size_t second_index)
      : quadratures_(&quadratures)
      , first_index_(first_index)
      , second_index_(second_index)
    {
      assert(first_index_ <= quadratures_->size());
      assert(second_index_ <= (*quadratures_)[first_index_].size());
    }

    MergedQuadratureIterator& operator++()
    {
      if (second_index_ != (*quadratures_)[first_index_].size() - 1) {
        ++second_index_;
      } else {
        // increase first_index_ until we reach either the next non-empty quadrature or the end
        while (++first_index_ < quadratures_->size() && !(*quadratures_)[first_index_].size())
          ;
        second_index_ = 0;
      }
      return *this;
    }

    MergedQuadratureIterator operator++(int)
    {
      auto ret = *this;
      ++(*this);
      return ret;
    }

    bool operator==(const MergedQuadratureIterator& other)
    {
      return (quadratures_ == other.quadratures_ && first_index_ == other.first_index_
              && second_index_ == other.second_index_);
    }

    bool operator!=(const MergedQuadratureIterator& other)
    {
      return !(*this == other);
    }

    const QuadPointType& operator*()
    {
      return (*quadratures_)[first_index_][second_index_];
    }

    size_t first_index()
    {
      return first_index_;
    }

    size_t second_index()
    {
      return first_index_;
    }

  private:
    const QuadraturesWrapperType* quadratures_;
    size_t first_index_;
    size_t second_index_;
  };

  using ConstIteratorType = MergedQuadratureIterator;

  ConstIteratorType begin() const
  {
    return ConstIteratorType(quadratures_, 0, 0);
  }

  ConstIteratorType end() const
  {
    return ConstIteratorType(quadratures_, quadratures_.size(), 0);
  }

  // iterator pointing to element at position ii in merged quadrature
  ConstIteratorType iterator(const size_t index) const
  {
    const auto indices = get_indices(index);
    return ConstIteratorType(quadratures_, indices[0], indices[1]);
  }

  const QuadPointType& operator[](const size_t index) const
  {
    const auto indices = get_indices(index);
    return quadratures_[indices[0]][indices[1]];
  }

  size_t size() const
  {
    return std::accumulate(
        quadratures_.begin(), quadratures_.end(), 0, [](const size_t& curr_size, const QuadratureType& quadrature) {
          return curr_size + quadrature.size();
        });
  }

private:
  FieldVector<size_t, 2> get_indices(size_t index) const
  {
    FieldVector<size_t, 2> ret{0, index};
    while (ret[1] >= quadratures_[ret[0]].size())
      ret[1] -= quadratures_[ret[0]++].size();
    return ret;
  }
  const QuadraturesWrapperType& quadratures_;
};

// Class to wrap several quadratures which can be used  either on their domain or as a single
// quadrature for the union of the individual domains.
template <class FieldType, size_t dimDomain>
class QuadraturesWrapper
{
public:
  using QuadratureType = Dune::QuadratureRule<FieldType, dimDomain>;
  using QuadPointType = Dune::QuadraturePoint<FieldType, dimDomain>;
  using QuadraturesType = std::vector<QuadratureType>;
  using IteratorType = typename QuadraturesType::iterator;
  using ConstIteratorType = typename QuadraturesType::const_iterator;
  using MergedQuadratureType = MergedQuadrature<FieldType, dimDomain>;

  QuadraturesWrapper(const size_t num_quadratures = 0, const QuadratureType& quadrature = QuadratureType())
    : quadratures_(num_quadratures, quadrature)
  {
  }

  QuadratureType& operator[](const size_t ii)
  {
    return quadratures_[ii];
  }

  const QuadratureType& operator[](const size_t ii) const
  {
    return quadratures_[ii];
  }

  IteratorType begin()
  {
    return quadratures_.begin();
  }

  IteratorType end()
  {
    return quadratures_.end();
  }

  ConstIteratorType begin() const
  {
    return quadratures_.begin();
  }

  ConstIteratorType end() const
  {
    return quadratures_.end();
  }

  MergedQuadratureType merged() const
  {
    return MergedQuadratureType(*this);
  }

  void clear()
  {
    quadratures_.clear();
  }

  size_t size() const
  {
    return quadratures_.size();
  }

private:
  QuadraturesType quadratures_;
};


} // namespace Problems
} // namespace Hyperbolic
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_TEST_HYPERBOLIC_SPHERICALQUADRATURES_WRAPPER_HH
