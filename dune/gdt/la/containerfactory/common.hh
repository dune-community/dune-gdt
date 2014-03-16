// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_LA_CONTAINERFACTORY_COMMON_HH
#define DUNE_GDT_LA_CONTAINERFACTORY_COMMON_HH

#include <dune/stuff/common/exceptions.hh>
#include <dune/stuff/la/container/common.hh>

#include "../containerfactory.hh"

namespace Dune {
namespace GDT {
namespace LA {


template< class ScalarType >
class ContainerFactory< Stuff::LA::CommonDenseMatrix< ScalarType > >
{
public:
  typedef Stuff::LA::CommonDenseMatrix< ScalarType > ContainerType;

  template< class T, class A >
  static ContainerType create(const SpaceInterface< T >& test_space, const SpaceInterface< A >& ansatz_space)
  {
    return ContainerType(test_space.mapper().size(), ansatz_space.mapper().size());
  } // ... create(...)

  template< class S >
  static ContainerType create(const SpaceInterface< S >& space)
  {
    return create(space, space);
  }

}; // class ContainerFactory


template< class ScalarType >
class ContainerFactory< Stuff::LA::CommonDenseVector< ScalarType > >
{
public:
  typedef Stuff::LA::CommonDenseVector< ScalarType > ContainerType;

  template< class T, class A >
  static ContainerType create(const SpaceInterface< T >& test_space, const SpaceInterface< A >& ansatz_space)
  {
    if (test_space.mapper.size() != ansatz_space.mapper().size())
      DUNE_THROW_COLORFULLY(Stuff::Exceptions::shapes_do_not_match,
                            "You called create(test_space, ansatz_space) for a Stuff::LA::CommonDenseVector "
                            << "with spaces of different sizes (" << test_space.mapper.size() << " and "
                            << ansatz_space.mapper.size() << ", respectively) which does not make any sense at all!");
    create(test_space);
  } // ... create(...)

  template< class S >
  static ContainerType create(const SpaceInterface< S >& space)
  {
    return ContainerType(space.mapper().size());
  }

}; // class ContainerFactory


} // namespace LA
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_LA_CONTAINERFACTORY_COMMON_HH
