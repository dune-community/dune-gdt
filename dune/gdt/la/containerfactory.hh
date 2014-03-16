// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LA_CONTAINERFACTORY_HH
#define DUNE_GDT_LA_CONTAINERFACTORY_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/common/type_utils.hh>

#include <dune/gdt/space/interface.hh>

namespace Dune {
namespace GDT {
namespace LA {


template <class ContainerImp>
class ContainerFactory
{
  static_assert(AlwaysFalse<ContainerImp>::value, "This is the unspecialized version of LA::ContainerFactory< ... >. "
                                                  "Please include the correct header for your container "
                                                  "implementation!");

public:
  typedef ContainerImp ContainerType;

  template <class T, class A>
  static ContainerType create(const SpaceInterface<T>& /*test_space*/, const SpaceInterface<A>& /*ansatz_space*/)
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::ContainerFactory< ... >. "
                          "Please include the correct header for your container implementation '"
                              << Stuff::Common::Typename<ContainerType>::value()
                              << "'!");
  } // ... create(...)

  template <class S>
  static ContainerType create(const SpaceInterface<S>& /*space*/)
  {
    DUNE_THROW_COLORFULLY(NotImplemented,
                          "This is the unspecialized version of LA::ContainerFactory< ... >. "
                          "Please include the correct header for your container implementation '"
                              << Stuff::Common::Typename<ContainerType>::value()
                              << "'!");
  } // ... create(...)

}; // class ContainerFactory


} // namespace LA
} // namespace GDT
} // namespace Dune

#include "containerfactory/common.hh"

#endif // DUNE_GDT_LA_CONTAINERFACTORY_HH
