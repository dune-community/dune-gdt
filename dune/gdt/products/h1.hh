// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_H1_HH
#define DUNE_GDT_PRODUCTS_H1_HH

#include "h1-internal.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief A localizable semi H1 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa LocalizableBase first and then \sa
 *        internal::H1SemiBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GridView, class Range, class Source = Range, class FieldType = double>
class H1SemiLocalizable : public LocalizableBase<internal::H1SemiBase<GridView, FieldType>, Range, Source>
{
  typedef LocalizableBase<internal::H1SemiBase<GridView, FieldType>, Range, Source> BaseType;

public:
  template <class... Args>
  H1SemiLocalizable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An assemblable semi H1 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::H1SemiBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class Matrix, class RangeSpace, class GridView = typename RangeSpace::GridViewType,
          class SourceSpace = RangeSpace, class FieldType = double>
class H1SemiAssemblable
    : public AssemblableBase<internal::H1SemiBase<GridView, FieldType>, Matrix, RangeSpace, SourceSpace>
{
  typedef AssemblableBase<internal::H1SemiBase<GridView, FieldType>, Matrix, RangeSpace, SourceSpace> BaseType;

public:
  template <class... Args>
  H1SemiAssemblable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief A semi H1 product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::H1SemiBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GridView, class FieldType = double>
class H1Semi : public GenericBase<internal::H1SemiBase<GridView, FieldType>>
{
  typedef GenericBase<internal::H1SemiBase<GridView, FieldType>> BaseType;

public:
  template <class... Args>
  H1Semi(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_H1_HH
