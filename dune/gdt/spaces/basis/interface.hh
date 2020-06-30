// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2013 - 2018)
//   Kirsten Weber   (2013)
//   René Fritze     (2014, 2016, 2018)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BASIS_INTERFACE_HH
#define DUNE_GDT_SPACES_BASIS_INTERFACE_HH

#include <dune/xt/grid/type_traits.hh>
#include <dune/xt/functions/interfaces/element-functions.hh>

#include <dune/gdt/local/finite-elements/interfaces.hh>

namespace Dune {
namespace GDT {


template <class E, size_t r = 1, size_t rC = 1, class R = double>
class LocalizedGlobalFiniteElementInterface : public XT::Functions::ElementFunctionSetInterface<E, r, rC, R>
{
  using BaseType = XT::Functions::ElementFunctionSetInterface<E, r, rC, R>;

public:
  using BaseType::d;
  using typename BaseType::D;
  using typename BaseType::DomainType;
  using typename BaseType::ElementType;
  using typename BaseType::RangeType;

  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC>;

  virtual ~LocalizedGlobalFiniteElementInterface() = default;

  /// \name ``These methods have to be implemented.''
  /// \{

  virtual const LocalFiniteElementType& finite_element() const = 0;

  virtual DynamicVector<R> default_data(const GeometryType& geometry_type) const = 0;

  virtual DynamicVector<R> backup() const = 0;

  virtual void restore(const ElementType& element, const DynamicVector<R>& data) = 0;

  virtual void interpolate(const std::function<RangeType(DomainType)>& element_function,
                           const int order,
                           DynamicVector<R>& dofs) const = 0;

  /// \}
  /// \name ``These methods are provided for convenience.''
  /// \{

  virtual void interpolate(const XT::Functions::ElementFunctionInterface<E, r, rC, R>& element_function,
                           DynamicVector<R>& dofs) const
  {
    this->interpolate([&](const auto& x) { return element_function.evaluate(x); }, element_function.order(), dofs);
  }

  template <class V, class GV>
  void interpolate(const std::function<RangeType(DomainType)>& element_function,
                   const int order,
                   LocalDofVector<V, GV>& dofs) const
  {
    const size_t sz = this->size();
    if (dofs_.size() < sz)
      dofs_.resize(sz);
    this->interpolate(element_function, order, dofs_);
    assert(dofs.size() >= sz);
    for (size_t ii = 0; ii < sz; ++ii)
      dofs.set_entry(ii, dofs_[ii]);
  } // ... interpolate(...)

  template <class V, class GV>
  void interpolate(const XT::Functions::ElementFunctionInterface<E, r, rC, R>& element_function,
                   LocalDofVector<V, GV>& dofs) const
  {
    const size_t sz = this->size();
    if (dofs_.size() < sz)
      dofs_.resize(sz);
    this->interpolate(element_function, dofs_);
    assert(dofs.size() >= sz);
    for (size_t ii = 0; ii < sz; ++ii)
      dofs.set_entry(ii, dofs_[ii]);
  } // ... interpolate(...)

  /// \}
  /// \name ``These methods are provided for convenience and should not be used within library code.''
  /// \{

  virtual DynamicVector<R> interpolate(const std::function<RangeType(DomainType)>& local_function,
                                       const int order) const
  {
    DynamicVector<R> result(this->size());
    this->interpolate(local_function, order, result);
    return result;
  }

  /// \}
private:
  mutable DynamicVector<R> dofs_;
}; // class LocalizedGlobalFiniteElementInterface


template <class GV, size_t range_dim = 1, size_t range_dim_columns = 1, class RangeField = double>
class GlobalBasisInterface
{
  static_assert(XT::Grid::is_view<GV>::value, "");

public:
  using GridViewType = GV;
  using E = XT::Grid::extract_entity_t<GridViewType>;
  using D = typename E::Geometry::ctype;
  static const constexpr size_t d = E::dimension;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using ElementType = E;
  using LocalizedType = LocalizedGlobalFiniteElementInterface<E, r, rC, R>;

  virtual ~GlobalBasisInterface() = default;

  virtual size_t max_size() const
  {
    return this->localize()->max_size();
  }

  virtual std::unique_ptr<LocalizedType> localize() const = 0;

  virtual std::unique_ptr<LocalizedType> localize(const ElementType& grid_element) const
  {
    auto fe = localize();
    fe->bind(grid_element);
    return fe;
  }

  /// \name These methods are required for grid adaptation.
  /// \{

  virtual void update_after_adapt()
  {
    DUNE_THROW(Exceptions::basis_error, "This basis does not support adaptation!");
  }

  /// \}
}; // class GlobalBasisInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASIS_INTERFACE_HH
