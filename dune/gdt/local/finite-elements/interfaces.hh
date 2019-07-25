// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)
//   René Fritze     (2018)
//   Tobias Leibner  (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH

#include <algorithm>
#include <functional>
#include <set>
#include <vector>

#include <dune/common/dynvector.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>

#include <dune/xt/common/unused.hh>
#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


// forward, include is below
template <class Vector, class GridView>
class LocalDofVector;


template <class DomainField, size_t domain_dim, class RangeField, size_t range_dim, size_t range_dim_columns = 1>
class LocalFiniteElementBasisInterface
{
  using ThisType = LocalFiniteElementBasisInterface<DomainField, domain_dim, RangeField, range_dim, range_dim_columns>;

public:
  using D = DomainField;
  static const constexpr size_t d = domain_dim;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using DomainType = FieldVector<D, d>;
  using RangeType = typename XT::Functions::RangeTypeSelector<R, r, rC>::type;
  using DerivativeRangeType = typename XT::Functions::DerivativeRangeTypeSelector<d, R, r, rC>::type;

  virtual ~LocalFiniteElementBasisInterface() = default;

  virtual ThisType* copy() const = 0;

  virtual const GeometryType& geometry_type() const = 0;

  virtual int order() const = 0;

  virtual size_t size() const = 0;

  virtual void evaluate(const DomainType& /*point_in_reference_element*/, std::vector<RangeType>& /*result*/) const = 0;

  virtual void jacobian(const DomainType& /*point_in_reference_element*/,
                        std::vector<DerivativeRangeType>& /*result*/) const = 0;

  /// \name ``These methods are provided for convenience and should not be used within library code.''
  /// \{

  virtual std::vector<RangeType> evaluate(const DomainType& point_in_reference_element) const
  {
    std::vector<RangeType> result(this->size());
    this->evaluate(point_in_reference_element, result);
    return result;
  }

  virtual std::vector<DerivativeRangeType> jacobian(const DomainType& point_in_reference_element) const
  {
    std::vector<DerivativeRangeType> result(this->size());
    this->jacobian(point_in_reference_element, result);
    return result;
  }

  /// \}
}; // class LocalFiniteElementBasisInterface


template <class DomainField, size_t domain_dim, class RangeField, size_t range_dim, size_t range_dim_columns = 1>
class LocalFiniteElementInterpolationInterface
{
  using ThisType =
      LocalFiniteElementInterpolationInterface<DomainField, domain_dim, RangeField, range_dim, range_dim_columns>;

public:
  using D = DomainField;
  static const constexpr size_t d = domain_dim;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using DomainType = FieldVector<D, d>;
  using RangeType = typename XT::Functions::RangeTypeSelector<R, r, rC>::type;

  LocalFiniteElementInterpolationInterface() {}

  virtual ~LocalFiniteElementInterpolationInterface() = default;

  virtual ThisType* copy() const = 0;

  virtual const GeometryType& geometry_type() const = 0;

  virtual size_t size() const = 0;

  virtual void interpolate(const std::function<RangeType(DomainType)>& local_function,
                           const int order,
                           DynamicVector<R>& dofs) const = 0;

  template <class V, class GV>
  void interpolate(const std::function<RangeType(DomainType)>& local_function,
                   const int order,
                   LocalDofVector<V, GV>& dofs) const
  {
    const size_t sz = this->size();
    if (dofs_.size() < sz)
      dofs_.resize(sz);
    this->interpolate(local_function, order, dofs_);
    for (size_t ii = 0; ii < sz; ++ii)
      dofs[ii] = dofs_[ii];
  } // ... interpolate(...)

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
}; // class LocalFiniteElementInterpolationInterface


template <class DomainField, size_t domain_dim>
class LocalFiniteElementCoefficientsInterface
{
  using ThisType = LocalFiniteElementCoefficientsInterface;

public:
  using D = DomainField;
  static const constexpr size_t d = domain_dim;

  virtual ~LocalFiniteElementCoefficientsInterface() = default;

  virtual ThisType* copy() const = 0;

  virtual const GeometryType& geometry_type() const = 0;

  virtual size_t size() const = 0;

  virtual const LocalKey& local_key(const size_t ii) const = 0;

  /**
   * \brief Computes the reverse information contained in the local keys.
   *
   * Returns a data structure data, such that data[codim][subentity_index] contains a (possibly empty) set of indices,
   * which can be used to access the corresponding local key:
\code
const auto codim_to_subentity_index_to_key_indices_map = local_key_indices(reference_element);
for (size_t codim = 0; codim < codim_to_subentity_index_to_key_indices_map.size(); ++codim) {
  const auto& subentity_index_to_key_indices_map = codim_to_subentity_index_to_key_indices_map[codim];
  for (size_t subentity_index = 0; subentity_index < subentity_index_to_key_indices_map.size(); ++subentity_index) {
    const auto& key_indices = subentity_index_to_key_indices_map[subentity_index];
    if (key_indices.size() == 0)
      std::cout << "no LocalKey associated with codim " << codim << " subentity " << subentity_index << std::endl;
    else {
      std::cout << "the following LocalKeys are associated with codim " << codim << " subentity " << subentity_index
                << ": " << std::endl;
      for (const auto& local_key_index : key_indices)
        std::cout << "  " << local_key(local_key_index) << std::endl;
    }
  }
}
\endcode
   *
   * \note It is guaranteed that access to map[codim][subentity_index] is valid for all 0 <= codim <= d and
   *       all 0 <= subentity_index < reference_element.size(codim).
   */
  std::vector<std::vector<std::vector<size_t>>> local_key_indices() const
  {
    // prepare
    const auto& reference_element = ReferenceElements<D, d>::general(geometry_type());
    std::vector<std::vector<std::vector<size_t>>> codim_to_subentity_index_to_key_indices_map(d + 1);
    for (int codim = 0; codim <= static_cast<int>(d); ++codim)
      codim_to_subentity_index_to_key_indices_map[codim] =
          std::vector<std::vector<size_t>>(reference_element.size(codim));
    // fill
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto& key = local_key(ii);
      auto& subentity_index_to_key_indices_map = codim_to_subentity_index_to_key_indices_map[key.codim()];
      subentity_index_to_key_indices_map[key.subEntity()].push_back(ii);
    }
    // sort
    for (auto& subentity_index_to_key_indices_map : codim_to_subentity_index_to_key_indices_map)
      for (auto& key_indices : subentity_index_to_key_indices_map) {
        std::sort(key_indices.begin(), key_indices.end());
        for (auto&& index : key_indices)
          DUNE_THROW_IF(
              std::count(key_indices.begin(), key_indices.end(), index) != 1, Exceptions::finite_element_error, "");
      }
    return codim_to_subentity_index_to_key_indices_map;
  } // ... local_key_indices(...)

  /**
   * \sa local_key_indices
   */
  std::vector<std::vector<size_t>> local_key_indices(const size_t codim) const
  {
    if (codim > d)
      DUNE_THROW(Exceptions::finite_element_error,
                 "d = " << d << "\n"
                        << "codim = " << codim);
    // pepare
    const auto& reference_element = ReferenceElements<D, d>::general(geometry_type());
    std::vector<std::vector<size_t>> subentity_index_to_key_indices_map(
        reference_element.size(XT::Common::numeric_cast<int>(codim)));
    // fill
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto& key = local_key(ii);
      if (key.codim() == codim)
        subentity_index_to_key_indices_map[key.subEntity()].push_back(ii);
    }
    // sort
    for (auto& key_indices : subentity_index_to_key_indices_map) {
      std::sort(key_indices.begin(), key_indices.end());
      for (auto&& index : key_indices)
        DUNE_THROW_IF(
            std::count(key_indices.begin(), key_indices.end(), index) != 1, Exceptions::finite_element_error, "");
    }
    return subentity_index_to_key_indices_map;
  } // ... local_key_indices(...)
}; // class LocalFiniteElementCoefficientsInterface


template <class DomainField, size_t domain_dim, class RangeField, size_t range_dim, size_t range_dim_columns = 1>
class LocalFiniteElementInterface
{
public:
  using D = DomainField;
  static const constexpr size_t d = domain_dim;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using DomainType = FieldVector<D, d>;
  using BasisType = LocalFiniteElementBasisInterface<D, d, R, r, rC>;
  using CoefficientsType = LocalFiniteElementCoefficientsInterface<D, d>;
  using InterpolationType = LocalFiniteElementInterpolationInterface<D, d, R, r, rC>;

  virtual ~LocalFiniteElementInterface() = default;

  virtual const GeometryType& geometry_type() const = 0;

  /**
   * \note This order does not have to coincide with basis().order()!
   *
   * E.g., for a local Raviart-Thomas finite element of order 0 on a simplex, the order of the finite element is 0, but
   * the order of the basis (a.k.a, its polynomial degree) is 1.
   */
  virtual int order() const = 0;

  virtual size_t size() const = 0;

  virtual const BasisType& basis() const = 0;

  virtual const CoefficientsType& coefficients() const = 0;

  virtual const InterpolationType& interpolation() const = 0;

  virtual bool is_lagrangian() const
  {
    return false;
  }

  /**
   * \brief Lagrange points associated with this local finite element (in reference element coordinates).
   *
   * Lagrange points are defined as usual, depending on the domain_dim and order() of the local finite element, and on
   * its range_dim or range_dim_cols.
   *
   * Thus it always holds that
   *
\code
assert(lagrange_points.size() == size() * r * rC);
\endcode
   *
   * Regarding the relation between a Lagrange point and an element of the basis, we have the following: in the scalar
   * case (r = rC = 1), it holds that
   *
\code
for (size_t ii = 0; ii < lagrange_points.size(); ++ii)
  for (size_t jj = 0; jj < basis().size(); ++jj)
    assert(basis().evaluate(lagrange_points()[ii])[jj] == (ii == jj ? 1 : 0));
\endcode
   *
   * while in the vector-valued case it holds that
   *
\code
for (size_t ii = 0; ii < lagrange_points.size(); ++ii)
  for (size_t rr = 0; rr < range_dim; ++rr)
    for (size_t jj = 0; jj < basis().size(); ++jj)
      assert(basis().evaluate(lagrange_points()[ii])[rr * lagrange_points.size() + jj][rr] == (ii == jj ? 1 : 0));
\endcode
   *
   * while the matrix-valued case is not implemented yet.
   *
   * See also LocalPowerFiniteElement for the reasoning behind the mapping in the vector-valued case.
   *
   * \sa LocalPowerFiniteElement
   */
  virtual const std::vector<DomainType>& lagrange_points() const
  {
    if (is_lagrangian())
      DUNE_THROW(NotImplemented,
                 "the implementor of this local finite element has to provider lagrange_points(), if "
                 "is_lagrangian() is true!");
    else
      DUNE_THROW(Exceptions::finite_element_error, "do not call lagrange_points() if is_lagrangian() is false!");
  }
}; // class LocalFiniteElementInterface


/// \todo allow for variable orders
template <class DomainField, size_t domain_dim, class RangeField, size_t range_dim, size_t range_dim_columns = 1>
class LocalFiniteElementFamilyInterface
{
public:
  using D = DomainField;
  static const constexpr size_t d = domain_dim;
  using R = RangeField;
  static const constexpr size_t r = range_dim;
  static const constexpr size_t rC = range_dim_columns;

  using LocalFiniteElementType = LocalFiniteElementInterface<D, d, R, r, rC>;

  virtual ~LocalFiniteElementFamilyInterface() = default;

  virtual const LocalFiniteElementType& get(const GeometryType& geometry_type, const int order) const = 0;
}; // class LocalFiniteElementFamilyInterface


} // namespace GDT
} // namespace Dune

#include <dune/gdt/local/dof-vector.hh>

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
