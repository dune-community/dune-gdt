// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)
//   Rene Milk       (2014)
//   Tobias Leibner  (2014)

#ifndef DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_FUNCTIONALS_INTERFACES_HH

#include <dune/common/deprecated.hh>

#include <dune/xt/common/crtp.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class FunctionalInterface : public Stuff::CRTPInterface<FunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::FieldType FieldType;

  template <class SourceType>
  FieldType apply(const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp().apply(source));
    return this->as_imp().apply(source);
  }
}; // class FunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_INTERFACES_HH
