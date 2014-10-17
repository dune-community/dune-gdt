// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_SWIPDGPENALTY_HH
#define DUNE_GDT_PRODUCTS_SWIPDGPENALTY_HH

#include "../../products/base.hh"
#include "swipdgpenalty-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {


/**
 * \brief A localizable product based on the SWIPDG penalty term.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa LocalizableBase first and then \sa
 *        internal::SwipdgPenaltyBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template< class GV, class FunctionImp, class TensorImp, class R, class S = R, class FieldType = double >
class SwipdgPenaltyLocalizable
  : public LocalizableBase< internal::SwipdgPenaltyBase< GV, FunctionImp, TensorImp, FieldType >, R, S >
{
  typedef LocalizableBase< internal::SwipdgPenaltyBase< GV, FunctionImp, TensorImp, FieldType >, R, S > BaseType;

public:
  template< class... Args >
  SwipdgPenaltyLocalizable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


/**
 * \brief An assemblable product based on the SWIPDG penalty term.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa AssemblableBase first and then \sa
 *        internal::SwipdgPenaltyBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template< class M, class F, class T, class R, class GV = typename R::GridViewType, class S = R, class FieldType = double >
class SwipdgPenaltyAssemblable
  : public AssemblableBase< internal::SwipdgPenaltyBase< GV, F, T, FieldType >, M, R, S >
{
  typedef AssemblableBase< internal::SwipdgPenaltyBase< GV, F, T, FieldType >, M, R, S > BaseType;

public:
  template< class... Args >
  SwipdgPenaltyAssemblable(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


/**
 * \brief A product based on the SWIPDG penalty term.
 *
 *        Possible ctor signaturer are a combination of the ones from \sa GenericBase first and then \sa
 *        internal::SwipdgPenaltyBase.
 * \todo  Add more documentation, especially a mathematical definition.
 */
template< class GV, class F, class T, class FieldType = double >
class SwipdgPenalty
  : public GenericBase< internal::SwipdgPenaltyBase< GV, F, T, FieldType > >
{
  typedef GenericBase< internal::SwipdgPenaltyBase< GV, F, T, FieldType > > BaseType;

public:
  template< class... Args >
  SwipdgPenalty(Args&& ...args)
    : BaseType(std::forward< Args >(args)...)
  {}
};


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_SWIPDGPENALTY_HH
