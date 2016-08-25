// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2013 - 2016)
//   Kirsten Weber   (2013)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_SPACES_BASEFUNCTIONSET_INTERFACE_HH
#define DUNE_GDT_SPACES_BASEFUNCTIONSET_INTERFACE_HH

#include <dune/xt/functions/interfaces.hh>
#include <dune/xt/common/crtp.hh>

namespace Dune {
namespace GDT {


/**
 *  \brief  The purpose of this interface is just to be used for template matching and to allow for access to the
 *          backend. All other functionality is enforced by XT::Functions::LocalfunctionSetInterface.
 *
 *          \see XT::Functions::LocalfunctionSetInterface for the template parameters D, d, R, r and rC.
 */
template <class Traits, class D, size_t d, class R, size_t r, size_t rC = 1>
class BaseFunctionSetInterface
    : public XT::Functions::LocalfunctionSetInterface<typename Traits::EntityType, D, d, R, r, rC>,
      public XT::CRTPInterface<BaseFunctionSetInterface<Traits, D, d, R, r, rC>, Traits>
{
  typedef XT::Functions::LocalfunctionSetInterface<typename Traits::EntityType, D, d, R, r, rC> BaseType;

public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::BackendType BackendType;
  typedef typename Traits::EntityType EntityType;

  explicit BaseFunctionSetInterface(const EntityType& ent)
    : BaseType(ent)
  {
  }

  const BackendType& backend() const
  {
    CHECK_CRTP(this->as_imp(*this).backend());
    return this->as_imp(*this).backend();
  }
}; // class BaseFunctionSetInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_BASEFUNCTIONSET_INTERFACE_HH
