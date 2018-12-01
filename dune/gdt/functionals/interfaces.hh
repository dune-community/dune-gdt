// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016, 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_FUNCTIONALS_INTERFACES_HH

#include <dune/xt/la/container/vector-interface.hh>

#include <dune/gdt/discretefunction/default.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


/**
 * \brief Interface for functionals.
 *
 * Consider a partition \tau_h (modelled by SourceGridView) of a physical domain, a (possibly vector- or matrix-values)
 * vector space K^{r \times rC} (modelled by Field, source_dim and source_dim_cols), this interface modells functionals
 *
 *   f: V_h -> K
 *
 * where V_h is the (discrete) space
 *
 *   V_h := \{ v_h: \tau_h -> K^{r \times rC} \}
 *
 * modelled by a SpaceInterface of appropriate dimensions. The functions v_h \in V_h, to which the funtional can be
 * applied, are modelled by ConstDiscreteFunction with an appropriate vector type derived from XT::LA::VectorInterface
 * (modelled by SourceVector).
 */
template <class SourceVector,
          class SourceGridView,
          size_t source_dim = 1,
          size_t source_dim_cols = 1,
          class Field = double>
class FunctionalInterface
{
public:
  static const constexpr size_t r = source_dim;
  static const constexpr size_t rC = source_dim_cols;
  using F = Field;

  using SourceSpaceType = SpaceInterface<SourceGridView, r, rC, F>;
  using SourceVectorType = SourceVector;
  using SourceType = ConstDiscreteFunction<SourceVectorType, SourceGridView, r, rC, F>;
  using FieldType = Field;

  virtual ~FunctionalInterface() = default;

  /// \name These methods have to be implemented.
  /// \{

  virtual bool linear() const = 0;

  virtual const SourceSpaceType& source_space() const = 0;

  virtual FieldType apply(const SourceType& source) const = 0;

  /// \}

  /**
   * Allows the implementation to do preparatory work (i.e., assemble the vector of a vector-based linear functional).
   *
   * \note In general, you have to call this method before calling apply!
   */
  virtual void assemble(const bool /*use_tbb*/ = 0) {}
}; // class FunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_INTERFACES_HH
