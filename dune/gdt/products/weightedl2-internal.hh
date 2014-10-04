// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_WEIGHTEDL2_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_WEIGHTEDL2_INTERNAL_HH

#include <type_traits>

#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localevaluation/product.hh>

#include "base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all weighted L2 products.
 * \note  Most likely you do not want to use this class directly, but Products::WeightedL2Localizable,
 *        Products::WeightedL2Assemblable or Products::WeightedL2 instead!
 */
template< class GV, class FunctionImp, class FieldImp >
class WeightedL2Base
  : public LocalOperatorProviderBase< GV >
{
  static_assert(Stuff::is_localizable_function< FunctionImp >::value,
                "FunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef WeightedL2Base< GV, FunctionImp, FieldImp > ThisType;
protected:
  typedef FunctionImp FunctionType;
public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral< LocalEvaluation::Product< FunctionType > > VolumeOperatorType;

  static const bool has_volume_operator = true;

  WeightedL2Base(const FunctionType& function, const size_t over_integrate = 0)
    : function_(function)
    , volume_operator_(over_integrate, function_)
  {}

  WeightedL2Base(const ThisType& other) = default;

private:
  const FunctionType& function_;
public:
  const VolumeOperatorType volume_operator_;
}; // WeightedL2Base


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_WEIGHTEDL2_INTERNAL_HH
