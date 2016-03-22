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
template< class GridView,
          class DiffusionFactor,
          class Range,
          class Source = Range,
          class FieldType = double,
          class DiffusionTensor = void >
class EllipticLocalizable
  : public LocalizableBase< internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor >,
                            Range, Source >
{
  typedef LocalizableBase
      < internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor >, Range, Source > BaseType;

public:
  template< class... Args >
  EllipticLocalizable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


/**
 * \brief An assemblable elliptic product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::EllipticBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template< class Matrix,
          class DiffusionFactor,
          class RangeSpace,
          class GridView = typename RangeSpace::GridViewType,
          class SourceSpace = RangeSpace,
          class FieldType = double,
          class DiffusionTensor = void >
class EllipticAssemblable
  : public AssemblableBase< internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor >,
                            Matrix, RangeSpace, SourceSpace >
{
  typedef AssemblableBase< internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor >,
                           Matrix, RangeSpace, SourceSpace > BaseType;

public:
  template< class... Args >
  EllipticAssemblable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


/**
 * \brief An elliptic product.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::EllipticBase. For instance:\code
Elliptic(grid_view, diffusion, over_integrate = 0);
Elliptic(grid_view, diffusion_factor, diffusion_tensor, over_integrate = 0);
\endcode
 * \todo  Add more documentation, especially a mathematical definition.
 */
template< class GridView, class DiffusionFactor, class FieldType = double, class DiffusionTensor = void >
class Elliptic
  : public GenericBase< internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor > >
{
  typedef GenericBase< internal::EllipticBase< DiffusionFactor, GridView, FieldType, DiffusionTensor > > BaseType;

public:
  template< class... Args >
  Elliptic(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


/**
 * \brief Creates a generic elliptic product.
 * \sa    Products::Elliptic
 *
 * Can be used, for instance, as\code
auto elliptic_product = Products::make_elliptic< double >(grid_view, diffusion_factor, diffusion_tensor, 1);
\endcode
 */
template< class RangeFieldType, class GV, class DF, class DT >
Elliptic< GV, DF, RangeFieldType, DT > make_elliptic(const GV& grid_view,
                                                     const DF& diffusion_factor,
                                                     const DT& diffusion_tensor, const size_t over_integrate = 0)
{
  return Elliptic< GV, DF, RangeFieldType, DT >(grid_view, diffusion_factor, diffusion_tensor, over_integrate);
}

/**
 * \brief Creates a generic elliptic product (variant with RangeFieldType = double).
 * \sa    make_elliptic
 * \sa    Products::Elliptic
 *
 * Can be used, for instance, as\code
auto elliptic_product = Products::make_elliptic(grid_view, diffusion_factor, diffusion_tensor, 1);
\endcode
 */
template< class GV, class DF, class DT >
Elliptic< GV, DF, double, DT > make_elliptic(const GV& grid_view,
                                             const DF& diffusion_factor,
                                             const DT& diffusion_tensor,
                                             const size_t over_integrate = 0)
{
  return make_elliptic< double >(grid_view, diffusion_factor, diffusion_tensor, over_integrate);
}


template< class GV, class R, class DF, class DT >
    EllipticLocalizable< GV, DF, R, R, typename DF::RangeFieldType, DT >
make_elliptic_localizable(const GV& grd_vw,
                          const R& rng,
                          const DF& diffusion_factor,
                          const DT& diffusion_tensor,
                          const size_t over_integrate = 0)
{
  typedef EllipticLocalizable< GV, DF, R, R, typename DF::RangeFieldType, DT > ProductType;
  return ProductType(grd_vw, rng, diffusion_factor, diffusion_tensor, over_integrate);
}


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ELLIPTIC_HH
