// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_H1_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_H1_INTERNAL_HH

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>

#include "elliptic-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all semi H1 products.
 *
 *        This class is implemented by using the seimplified variant of EllipticBase with a Stuff::Functions::Constant
 *        of value 1 as the diffusion.
 * \note  Most likely you do not want to use this class directly, but Products::H1SemiLocalizable,
 *        Products::H1SemiAssemblable or Products::H1Semi instead!
 */
template< class GV, class FieldImp >
class H1SemiBase
  : DSC::ConstStorageProvider
    < Stuff::Functions::Constant
          < typename GV::template Codim< 0 >::Entity, typename GV::ctype, GV::dimension, FieldImp, 1 > >
  , public EllipticBase< Stuff::Functions::Constant
          < typename GV::template Codim< 0 >::Entity, typename GV::ctype, GV::dimension, FieldImp, 1 >, GV, FieldImp >
{
  typedef DSC::ConstStorageProvider< Stuff::Functions::Constant
      < typename GV::template Codim< 0 >::Entity, typename GV::ctype, GV::dimension, FieldImp, 1 > >
                                     StorageBaseType;
  typedef EllipticBase< Stuff::Functions::Constant
      < typename GV::template Codim< 0 >::Entity, typename GV::ctype, GV::dimension, FieldImp, 1 >, GV, FieldImp >
                                     EllipticBaseType;
  typedef H1SemiBase< GV, FieldImp > ThisType;
public:

  H1SemiBase(const size_t over_integrate = 0)
    : StorageBaseType(new typename EllipticBaseType::DiffusionType(1))
    , EllipticBaseType(this->storage_access(), over_integrate)
    , over_integrate_(over_integrate)
  {}

  /**
   * \note We need the manual copy ctor bc of the Stuff::Common::ConstStorageProvider
   */
  H1SemiBase(const ThisType& other)
    : StorageBaseType(new typename EllipticBaseType::DiffusionType(1))
    , EllipticBaseType(this->storage_access(), other.over_integrate_)
    , over_integrate_(other.over_integrate_)
  {}

private:
  const size_t over_integrate_; //!< needed to provide manual copy ctor
}; // class H1SemiBase


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_H1_INTERNAL_HH
