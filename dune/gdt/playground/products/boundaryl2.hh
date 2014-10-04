// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BOUNDARYL2_HH
#define DUNE_GDT_PRODUCTS_BOUNDARYL2_HH

#include "../../products/base.hh"
#include "boundaryl2-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief A localizable L2 product over the boundary of a domain.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa LocalizableBase first and then \sa
 *        internal::BoundaryL2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GV, class R, class S = R, class FieldType = double>
class BoundaryL2Localizable : public LocalizableBase<internal::BoundaryL2Base<GV, FieldType>, R, S>
{
  typedef LocalizableBase<internal::BoundaryL2Base<GV, FieldType>, R, S> BaseType;

public:
  template <class... Args>
  BoundaryL2Localizable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An assemblable L2 product over the boundary of a domain.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::BoundaryL2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class M, class R, class GV = typename R::GridViewType, class S = R, class FieldType = double>
class BoundaryL2Assemblable : public AssemblableBase<internal::BoundaryL2Base<GV, FieldType>, M, R, S>
{
  typedef AssemblableBase<internal::BoundaryL2Base<GV, FieldType>, M, R, S> BaseType;

public:
  template <class... Args>
  BoundaryL2Assemblable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An L2 product over the boundary of a domain.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::L2Base.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GV, class FieldType = double>
class BoundaryL2 : public GenericBase<internal::BoundaryL2Base<GV, FieldType>>
{
  typedef GenericBase<internal::BoundaryL2Base<GV, FieldType>> BaseType;

public:
  template <class... Args>
  BoundaryL2(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BOUNDARYL2_HH
