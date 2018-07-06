// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2014, 2016, 2018)
//   Tobias Leibner  (2014, 2017)

#ifndef DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_FUNCTIONALS_INTERFACES_HH

#include <dune/xt/common/crtp.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class FunctionalInterface : public XT::Common::CRTPInterface<FunctionalInterface<Traits>, Traits>
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
