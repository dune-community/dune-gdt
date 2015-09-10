// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_HH

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/common/configuration.hh>
#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container.hh>

#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/integrals.hh>
#include <dune/gdt/operators/interfaces.hh>
#include <dune/gdt/operators/default.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {


// ////////////////////////// //
// EllipticLocalizableProduct //
// ////////////////////////// //

template< class DiffusionFactorType,
          typename DiffusionTensorType, // may be void
          class GridView,
          class Range,
          class Source = Range,
          class Field = typename Range::RangeFieldType >
class EllipticLocalizableProduct
  : public LocalizableProductDefault< GridView, Range, Source, Field >
{
  typedef LocalizableProductDefault< GridView, Range, Source, Field >           BaseType;
  typedef LocalVolumeIntegralOperator
      < LocalEvaluation::Elliptic< DiffusionFactorType, DiffusionTensorType > > LocalEllipticOperatorType;
public:
  // Usually, we only hace to hold the data functions for the local operator and perfect forward the rest of the
  // arguments to BaseType. Here, it is a bit more complicated, since DiffusionTensorType might be void (and
  // DiffusionFactorType is the the only diffusion). To handle this case we require the enable_if hacks below to
  // disable half of the ctors if DiffusionTensorType is void. In addition we require each ctor twice, once with
  // over_integrate, once without (since the perfect forwarding does not allow for default arguments).

  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value) //
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const DiffusionImp& diffusion, Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionImp // This ctor is only enabled if we are given a single diffusion data function.
          , typename = typename std::enable_if<    (std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const size_t over_integrate,
                                      const DiffusionImp& diffusion,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

  template< typename DiffusionFactorImp // This ctor is only enabled
          , typename DiffusionTensorImp // if we are given two diffusion data functions (factor and tensor).
          , typename = typename std::enable_if<    (!std::is_same< DiffusionTensorType, void >::value)
                                                && (std::is_same< DiffusionFactorImp, DiffusionFactorType >::value)
                                                && sizeof(DiffusionFactorImp) >::type
          , class ...Args >
  explicit EllipticLocalizableProduct(const size_t over_integrate,
                                      const DiffusionFactorImp& diffusion_factor,
                                      const DiffusionTensorImp& diffusion_tensor,
                                      Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
    , local_elliptic_operator_(over_integrate, diffusion_factor, diffusion_tensor)
  {
    this->add(local_elliptic_operator_);
  }

private:
  const LocalEllipticOperatorType local_elliptic_operator_;
}; // class EllipticLocalizableProduct


// ///////////////////////////////// //
// make_elliptic_localizable_product //
// ///////////////////////////////// //

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion.
 */
template< class DiffusionType, class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                            , std::unique_ptr< EllipticLocalizableProduct< DiffusionType, void, GridViewType,
                                                                           RangeType, SourceType > >
                            >::type
make_elliptic_localizable_product(const DiffusionType& diffusion,
                                  const GridViewType& grid_view,
                                  const RangeType& range,
                                  const SourceType& source,
                                  const size_t over_integrate = 0)
{
  return DSC::make_unique
      < EllipticLocalizableProduct< DiffusionType, void, GridViewType, RangeType, SourceType > >(
          over_integrate, diffusion, grid_view, range, source);
}

/**
 * \sa EllipticLocalizableProduct, especially for the role of diffusion_factor and diffusion_tensor.
 */
template< class DiffusionFactorType, class DiffusionTensorType, class GridViewType, class RangeType, class SourceType >
    typename std::enable_if<    Stuff::is_localizable_function< DiffusionFactorType >::value
                             && Stuff::is_localizable_function< DiffusionTensorType >::value
                             && Stuff::is_localizable_function< RangeType >::value
                             && Stuff::is_localizable_function< SourceType >::value
                           , std::unique_ptr< EllipticLocalizableProduct< DiffusionFactorType, DiffusionTensorType,
                                                                          GridViewType, RangeType, SourceType > >
                           >::type
make_elliptic_localizable_product(const DiffusionFactorType& diffusion_factor,
                                  const DiffusionTensorType& diffusion_tensor,
                                  const GridViewType& grid_view,
                                  const RangeType& range,
                                  const SourceType& source,
                                  const size_t over_integrate = 0)
{
  return DSC::make_unique
      < EllipticLocalizableProduct< DiffusionFactorType, DiffusionTensorType, GridViewType, RangeType, SourceType > >(
          over_integrate, diffusion_factor, diffusion_tensor, grid_view, range, source);
}


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_HH
