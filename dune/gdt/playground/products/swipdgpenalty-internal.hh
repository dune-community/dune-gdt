// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_SWIPDGPENALTY_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_SWIPDGPENALTY_INTERNAL_HH

#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/localoperator/codim1.hh>

#include "../../products/base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace internal {


/**
 * \brief Base class for all products based on the SWIPDG penalty term.
 * \note  Most likely you do not want to use this class directly, but Products::SwipdgPenaltyLocalizable,
 *        Products::SwipdgPenaltyAssemblable or Products::SwipdgPenalty instead!
 */
template< class GV, class FunctionImp, class TensorImp, class FieldImp >
class SwipdgPenaltyBase
  : public LocalOperatorProviderBase< GV >
{
  static_assert(Stuff::is_localizable_function< FunctionImp >::value,
                "FunctionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< TensorImp >::value,
                "TensorImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef SwipdgPenaltyBase< GV, FunctionImp, TensorImp, FieldImp > ThisType;
public:
  typedef GV          GridViewType;
  typedef FunctionImp FunctionType;
  typedef TensorImp   TensorType;
  typedef FieldImp    FieldType;

  typedef LocalOperator::Codim1CouplingIntegral
      < LocalEvaluation::SWIPDG::InnerPenalty< FunctionType, TensorType > >       CouplingOperatorType;
  typedef LocalOperator::Codim1BoundaryIntegral
      < LocalEvaluation::SWIPDG::BoundaryLHSPenalty< FunctionType, TensorType > > BoundaryOperatorType;

  static const bool has_coupling_operator = true;
  static const bool has_boundary_operator = true;

  SwipdgPenaltyBase(const FunctionType& func, const TensorType& tensor, const size_t over_integrate = 0)
    : coupling_operator_(over_integrate, func, tensor)
    , boundary_operator_(over_integrate, func, tensor)
  {}

  SwipdgPenaltyBase(const ThisType& other) = default;

  const CouplingOperatorType coupling_operator_;
  const BoundaryOperatorType boundary_operator_;
}; // SwipdgPenaltyBase


} // namespace internal
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_SWIPDGPENALTY_INTERNAL_HH
