// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH

#include <functional>
#include <set>
#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/localkey.hh>

#include <dune/xt/functions/type_traits.hh>

#include <dune/gdt/exceptions.hh>

namespace Dune {
namespace GDT {


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

  /**
   * \name ``These methods are provided for convenience and should not be used within library code.''
   */
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

  virtual ~LocalFiniteElementInterpolationInterface() = default;

  virtual ThisType* copy() const = 0;

  virtual const GeometryType& geometry_type() const = 0;

  virtual size_t size() const = 0;

  virtual void interpolate(const std::function<RangeType(DomainType)>& local_function,
                           const int order,
                           std::vector<R>& dofs) const = 0;

  /// \name ``These methods are provided for convenience and should not be used within library code.''
  /// \{

  virtual std::vector<R> interpolate(const std::function<RangeType(DomainType)>& local_function, const int order) const
  {
    std::vector<R> result(this->size());
    this->interpolate(local_function, order, result);
    return result;
  }

  /// \}
}; // class LocalFiniteElementInterpolationInterface


template <class DomainField, size_t domain_dim>
class LocalFiniteElementCoefficientsInterface
{
  using ThisType = LocalFiniteElementCoefficientsInterface<DomainField, domain_dim>;

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
  std::vector<std::vector<std::set<size_t>>> local_key_indices() const
  {
    const auto& reference_element = ReferenceElements<D, d>::general(geometry_type());
    // pepare
    std::vector<std::vector<std::set<size_t>>> codim_to_subentity_index_to_key_indices_map(d + 1);
    for (size_t codim = 0; codim <= d; ++codim)
      codim_to_subentity_index_to_key_indices_map[codim] = std::vector<std::set<size_t>>(reference_element.size(codim));
    // fill
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto& key = local_key(ii);
      auto& subentity_index_to_key_indices_map = codim_to_subentity_index_to_key_indices_map[key.codim()];
      subentity_index_to_key_indices_map[key.subEntity()].insert(ii);
    }
    return codim_to_subentity_index_to_key_indices_map;
  } // ... local_key_indices(...)

  /**
   * \sa local_key_indices
   */
  std::vector<std::set<size_t>> local_key_indices(const size_t codim) const
  {
    const auto& reference_element = ReferenceElements<D, d>::general(geometry_type());
    if (codim > d)
      DUNE_THROW(Exceptions::finite_element_error,
                 "d = " << d << "\n"
                        << "codim = "
                        << codim);
    // pepare
    std::vector<std::set<size_t>> subentity_index_to_key_indices_map(reference_element.size(codim));
    // fill
    for (size_t ii = 0; ii < size(); ++ii) {
      const auto& key = local_key(ii);
      if (key.codim() == codim)
        subentity_index_to_key_indices_map[key.subEntity()].insert(ii);
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


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_INTERFACES_HH
