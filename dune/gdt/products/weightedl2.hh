// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_WEIGHTEDL2_HH
#define DUNE_GDT_PRODUCTS_WEIGHTEDL2_HH

#include "base.hh"
#include "weightedl2-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief A localizable weighted L2 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa LocalizableBase first and then \sa
 *        internal::WeightedL2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GV, class F, class R, class S = R, class FieldType = double>
class WeightedL2Localizable : public LocalizableBase<internal::WeightedL2Base<GV, F, FieldType>, R, S>
{
  typedef LocalizableBase<internal::WeightedL2Base<GV, F, FieldType>, R, S> BaseType;

public:
  template <class... Args>
  WeightedL2Localizable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An assemblable weighted L2 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::WeightedL2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class M, class F, class R, class GV = typename R::GridViewType, class S = R, class FieldType = double>
class WeightedL2Assemblable : public AssemblableBase<internal::WeightedL2Base<GV, F, FieldType>, M, R, S>
{
  typedef AssemblableBase<internal::WeightedL2Base<GV, F, FieldType>, M, R, S> BaseType;

public:
  template <class... Args>
  WeightedL2Assemblable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief A weighted L2 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::WeightedL2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GV, class F, class FieldType = double>
class WeightedL2 : public GenericBase<internal::WeightedL2Base<GV, F, FieldType>>
{
  typedef GenericBase<internal::WeightedL2Base<GV, F, FieldType>> BaseType;

public:
  template <class... Args>
  WeightedL2(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_WEIGHTEDL2_HH
