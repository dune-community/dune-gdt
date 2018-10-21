// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2017)
//   Rene Milk       (2014, 2016 - 2018)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_SPACES_MAPPER_INTERFACES_HH
#define DUNE_GDT_SPACES_MAPPER_INTERFACES_HH

#include <dune/common/dynvector.hh>

#include <dune/xt/grid/type_traits.hh>

#include <dune/gdt/local/finite-elements/interfaces.hh>

namespace Dune {
namespace GDT {


template <class GV>
class MapperInterface
{
  static_assert(XT::Grid::is_view<GV>::value, "");

public:
  using GridViewType = GV;
  using D = typename GV::ctype;
  static const constexpr size_t d = GV::dimension;
  using ElementType = XT::Grid::extract_entity_t<GridViewType>;

  virtual ~MapperInterface() = default;

  virtual const GridViewType& grid_view() const = 0;

  virtual const LocalFiniteElementCoefficientsInterface<D, d>&
  local_coefficients(const GeometryType& geometry_type) const = 0;

  virtual size_t size() const = 0;

  virtual size_t max_local_size() const = 0;

  virtual size_t local_size(const ElementType& element) const = 0;

  virtual size_t global_index(const ElementType& element, const size_t local_index) const = 0;

  virtual void global_indices(const ElementType& element, DynamicVector<size_t>& indices) const = 0;

  virtual DynamicVector<size_t> global_indices(const ElementType& element) const
  {
    DynamicVector<size_t> ret(local_size(element));
    global_indices(element, ret);
    return ret;
  }

  /// \name These methods are required for grid adaptation.
  /// \{

  virtual void update_after_adapt()
  {
    DUNE_THROW(Exceptions::mapper_error, "This mapper does not support adaptation!");
  }

  /// \}
}; // class MapperInterface


} // namespace GDT
} // namespace Dune


#endif // DUNE_GDT_SPACES_MAPPER_INTERFACES_HH
