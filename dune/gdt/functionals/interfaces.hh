// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_FUNCTIONALS_INTERFACES_HH

#include <dune/common/deprecated.hh>

#include <dune/stuff/common/crtp.hh>

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
