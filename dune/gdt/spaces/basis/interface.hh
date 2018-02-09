// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Kirsten Weber   (2013)
//   Rene Milk       (2014, 2016)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BASIS_INTERFACE_HH
#define DUNE_GDT_SPACES_BASIS_INTERFACE_HH

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/local-functions.hh>

#include <dune/gdt/local/finite-elements/interfaces.hh>

namespace Dune {
namespace GDT {


template <class Element, size_t range_dim = 1, size_t range_dim_columns = 1, class RangeField = double>
class GlobalBasisInterface
{
  static_assert(XT::Grid::is_entity<Element>::value, "");

public:
  using E = Element;
  using D = typename E::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using ElementType = E;
  using ShapeFunctionsType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;
  using LocalizedBasisType = XT::Functions::LocalfunctionSetInterface<E, D, d, R, r, rC>;

  virtual ~GlobalBasisInterface() = default;

  virtual const ShapeFunctionsType& shape_functions(const GeometryType& geometry_type) const = 0;

  virtual std::unique_ptr<LocalizedBasisType> localize(const ElementType& element) const = 0;
}; // class GlobalBasisInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_INTERFACE_HH
