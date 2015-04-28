// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//
// Contributors: Tobias Leibner

#ifndef DUNE_GDT_SPACES_PRODUCTINTERFACE_HH
#define DUNE_GDT_SPACES_PRODUCTINTERFACE_HH

#include <tuple>

#include "interface.hh"

namespace Dune {
namespace GDT {

template< class Traits >
class ProductSpaceInterface
    : public SpaceInterface< Traits, Traits::dimDomain, Traits::dimRange, Traits::dimRangeCols >
{
  typedef typename Traits::FactorMapperType FactorMapperType;
  typedef typename Traits::SpaceTupleType SpaceTupleType;

public:
  /**
   * \defgroup interface ´´These methods have to be implemented in addition to the SpaceInterface methods!''
   * @{
   **/
  const FactorMapperType& factor_mapper() const
  {
    CHECK_CRTP(this->as_imp().factor_mapper());
    return this->as_imp().factor_mapper();
  }

  template< size_t i >
  const typename std::tuple_element< i, SpaceTupleType >::type& factor() const
  {
    CHECK_CRTP(this->as_imp().template factor< i >());
    return this->as_imp().template factor< i >();
  }
  /**
    }
   **/

}; // class ProductSpaceInterface< Traits >


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_SPACES_PRODUCTINTERFACE_HH
