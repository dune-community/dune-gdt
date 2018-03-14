// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2018)

#ifndef DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH
#define DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH

#include <dune/xt/common/memory.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {


/**
 * \brief Implements LocalFiniteElementInterface, given a basis, coefficients and an interpolation.
 */
template <class D, size_t d, class R, size_t r, size_t rC = 1>
class LocalFiniteElementDefault : public LocalFiniteElementInterface<D, d, R, r, rC>
{
  using ThisType = LocalFiniteElementDefault<D, d, R, r, rC>;
  using BaseType = LocalFiniteElementInterface<D, d, R, r, rC>;

public:
  using typename BaseType::DomainType;
  using typename BaseType::BasisType;
  using typename BaseType::CoefficientsType;
  using typename BaseType::InterpolationType;

  LocalFiniteElementDefault(const int ord,
                            const BasisType& bas,
                            const CoefficientsType& coeffs,
                            const InterpolationType& inter,
                            std::vector<DomainType> lps = {})
    : geometry_type_(bas.geometry_type())
    , order_(ord)
    , basis_(bas)
    , coefficients_(coeffs)
    , interpolation_(inter)
    , lagrange_points_(lps)
  {
    DUNE_THROW_IF(coeffs.geometry_type() != geometry_type_, Exceptions::finite_element_error, "");
    DUNE_THROW_IF(coeffs.size() != bas.size(), Exceptions::finite_element_error, "");
    DUNE_THROW_IF(inter.geometry_type() != geometry_type_, Exceptions::finite_element_error, "");
    DUNE_THROW_IF(inter.size() != bas.size(), Exceptions::finite_element_error, "");
    DUNE_THROW_IF(!(lps.size() == bas.size() || lps.size() == 0), Exceptions::finite_element_error, "");
  }

  /**
    * \attention Do not delete any of the moved in raw pointers manually afterwads!
    **/
  LocalFiniteElementDefault(const int ord,
                            BasisType*&& bas_ptr,
                            CoefficientsType*&& coeffs_ptr,
                            InterpolationType*&& inter_ptr,
                            std::vector<DomainType> lps = {})
    : geometry_type_(bas_ptr->geometry_type())
    , order_(ord)
    , basis_(std::move(bas_ptr))
    , coefficients_(std::move(coeffs_ptr))
    , interpolation_(std::move(inter_ptr))
    , lagrange_points_(lps)
  {
    DUNE_THROW_IF(coeffs_ptr->geometry_type() != geometry_type_, Exceptions::finite_element_error, "");
    DUNE_THROW_IF(coeffs_ptr->size() != bas_ptr->size(), Exceptions::finite_element_error, "");
    DUNE_THROW_IF(inter_ptr->geometry_type() != geometry_type_, Exceptions::finite_element_error, "");
    DUNE_THROW_IF(inter_ptr->size() != bas_ptr->size(), Exceptions::finite_element_error, "");
    DUNE_THROW_IF(!(lps.size() == bas_ptr->size() || lps.size() == 0), Exceptions::finite_element_error, "");
  }

  const GeometryType& geometry_type() const override
  {
    return geometry_type_;
  }

  int order() const override
  {
    return order_;
  }

  size_t size() const override
  {
    return basis_.access().size();
  }

  const BasisType& basis() const override
  {
    return basis_.access();
  }

  const CoefficientsType& coefficients() const override
  {
    return coefficients_.access();
  }

  const InterpolationType& interpolation() const override
  {
    return interpolation_.access();
  }

  bool is_lagrangian() const override
  {
    return lagrange_points_.size() > 0;
  }

  const std::vector<DomainType>& lagrange_points() const override
  {
    if (!is_lagrangian())
      DUNE_THROW(Exceptions::finite_element_error, "do not call lagrange_points() if is_lagrangian() is false!");
    return lagrange_points_;
  }

private:
  const GeometryType geometry_type_;
  const int order_;
  const XT::Common::ConstStorageProvider<BasisType> basis_;
  const XT::Common::ConstStorageProvider<CoefficientsType> coefficients_;
  const XT::Common::ConstStorageProvider<InterpolationType> interpolation_;
  const std::vector<DomainType> lagrange_points_;
}; // class LocalFiniteElementDefault


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCAL_FINITE_ELEMENTS_DEFAULT_HH
