// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ELLIPTIC_HH
#define DUNE_GDT_PRODUCTS_ELLIPTIC_HH

#include "elliptic-internal.hh"
#include "base.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief A localizable elliptic product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa LocalizableBase first and then \sa
 *        internal::EllipticBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GridView, class DiffusionFactor, class Range, class Source = Range, class FieldType = double,
          class DiffusionTensor = void>
class EllipticLocalizable
    : public LocalizableBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>, Range,
                             Source>
{
  typedef LocalizableBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>, Range, Source>
      BaseType;

public:
  template <class... Args>
  EllipticLocalizable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An assemblable elliptic product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::EllipticBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class Matrix, class DiffusionFactor, class RangeSpace, class GridView = typename RangeSpace::GridViewType,
          class SourceSpace = RangeSpace, class FieldType = double, class DiffusionTensor = void>
class EllipticAssemblable
    : public AssemblableBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>, Matrix,
                             RangeSpace, SourceSpace>
{
  typedef AssemblableBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>, Matrix,
                          RangeSpace, SourceSpace> BaseType;

public:
  template <class... Args>
  EllipticAssemblable(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


/**
 * \brief An elliptic product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::EllipticBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template <class GridView, class DiffusionFactor, class FieldType = double, class DiffusionTensor = void>
class Elliptic : public GenericBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>>
{
  typedef GenericBase<internal::EllipticBase<DiffusionFactor, GridView, FieldType, DiffusionTensor>> BaseType;

public:
  template <class... Args>
  Elliptic(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ELLIPTIC_HH
