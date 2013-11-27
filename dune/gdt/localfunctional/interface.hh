// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH
#define DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH

#include <vector>

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/common/dynvector.hh>

#include <dune/gdt/basefunctionset/interface.hh>

namespace Dune {
namespace GDT {
namespace LocalFunctional {


template <class Traits>
class Codim0Interface
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numTmpObjectsRequired());
    return asImp().numTmpObjectsRequired();
  }

  /**
   *  \brief      Applies the local functional.
   *  \tparam T   Traits of the test BaseFunctionSetInterface
   *  \tparam D   DomainFieldType
   *  \tparam d   dimDomain
   *  \tparam R   RangeFieldType
   *  \tparam r   dimRange
   *  \tparam rC  dimRangeCols
   */
  template <class T, class D, int d, class R, int r, int rC>
  void apply(const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase, Dune::DynamicVector<R>& ret,
             std::vector<Dune::DynamicVector<R>>& tmpLocalVectors) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().apply(testBase, ret, tmpLocalVectors));
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Codim0Interface

template <class Traits>
class Codim1Interface
{
public:
  typedef typename Traits::derived_type derived_type;

  size_t numTmpObjectsRequired() const
  {
    CHECK_INTERFACE_IMPLEMENTATION(asImp().numTmpObjectsRequired());
    return asImp().numTmpObjectsRequired();
  }

  /**
   *  \brief Applies the local functional.
   *  \tparam T                 Traits of the test BaseFunctionSetInterface implementation
   *  \tparam IntersectionType  A model of Dune::Intersection< ... >
   *  \tparam D                 DomainFieldType
   *  \tparam d                 dimDomain
   *  \tparam R                 RangeFieldType
   *  \tparam r                 dimRange of the of the testBase
   *  \tparam rC                dimRangeCols of the testBase
   */
  template <class T, class IntersectionType, class D, int d, class R, int r, int rC>
  void apply(const BaseFunctionSetInterface<T, D, d, R, r, rC>& testBase, const IntersectionType& intersection,
             Dune::DynamicVector<R>& ret, std::vector<Dune::DynamicVector<R>>& tmpLocalVectors) const
  {
    CHECK_AND_CALL_INTERFACE_IMPLEMENTATION(asImp().apply(testBase, intersection, ret, tmpLocalVectors));
  }

  derived_type& asImp()
  {
    return static_cast<derived_type&>(*this);
  }

  const derived_type& asImp() const
  {
    return static_cast<const derived_type&>(*this);
  }
}; // class Codim1Interface


} // namespace LocalFunctional
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LOCALFUNCTIONAL_INTERFACE_HH
