// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_BOUNDARYL2_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_BOUNDARYL2_INTERNAL_HH

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/constant.hh>

#include <dune/gdt/localevaluation/product.hh>
#include <dune/gdt/localoperator/codim1.hh>

#include "../../products/base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all L2 products over the boundary of a domain.
 *
 * \todo  The boundary integral is computed over all boundary intersections atm, but we might add an optional
 *        boundary info and integrate over dirichlet intersections only.
 */
template <class GV, class FieldImp>
class BoundaryL2Base
    : DSC::ConstStorageProvider<Stuff::Functions::Constant<typename GV::template Codim<0>::Entity, typename GV::ctype,
                                                           GV::dimension, FieldImp, 1>>,
      public LocalOperatorProviderBase<GV>
{
  typedef DSC::ConstStorageProvider<Stuff::Functions::Constant<
      typename GV::template Codim<0>::Entity, typename GV::ctype, GV::dimension, FieldImp, 1>> FunctionProvider;
  typedef BoundaryL2Base<GV, FieldImp> ThisType;
  typedef Stuff::Functions::Constant<typename GV::template Codim<0>::Entity, typename GV::ctype, GV::dimension,
                                     FieldImp, 1> FunctionType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;

  typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::Product<FunctionType>> BoundaryOperatorType;

  static const bool has_boundary_operator = true;

  BoundaryL2Base(const size_t over_integrate = 0)
    : FunctionProvider(new FunctionType(1))
    , over_integrate_(over_integrate)
    , boundary_operator_(over_integrate_, FunctionProvider::storage_access())
  {
  }

  /**
   * \note We need the manual copy ctor bc of the Stuff::Common::ConstStorageProvider
   */
  BoundaryL2Base(const ThisType& other)
    : FunctionProvider(new FunctionType(1))
    , over_integrate_(other.over_integrate_)
    , boundary_operator_(other.over_integrate_, FunctionProvider::storage_access())
  {
  }

private:
  const size_t over_integrate_; //!< needed to provide manual copy ctor
public:
  const BoundaryOperatorType boundary_operator_;
}; // BoundaryL2Base


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_BOUNDARYL2_INTERNAL_HH
