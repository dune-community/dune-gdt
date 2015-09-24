// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_INTERFACES_HH
#define DUNE_GDT_FUNCTIONALS_INTERFACES_HH

#include <dune/stuff/common/crtp.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/discretefunction/default.hh>

namespace Dune {
namespace GDT {


template< class Traits >
class FunctionalInterface
  : public Stuff::CRTPInterface< FunctionalInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::FieldType    FieldType;

  template< class SourceType >
  FieldType apply(const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp().apply(source));
    return this->as_imp().apply(source);
  }
}; // class FunctionalInterface


//! \note derive from FunctionalInterface
template< class Traits >
class AssemblableFunctionalInterface
  : protected Stuff::CRTPInterface< AssemblableFunctionalInterface< Traits >, Traits >
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SpaceType    SpaceType;
  typedef typename Traits::VectorType   VectorType;
  typedef typename Traits::ScalarType   ScalarType;

private:
  static_assert(is_space< SpaceType >::value, "SpaceType has to be derived from SpaceInterface!");
  static_assert(Stuff::LA::is_vector< VectorType >::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");

public:
  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp().grid_view());
    return this->as_imp().grid_view();
  }

  const SpaceType& space() const
  {
    CHECK_CRTP(this->as_imp().space());
    return this->as_imp().space();
  }

  void assemble()
  {
    CHECK_AND_CALL_CRTP(this->as_imp().assemble());
  }

  VectorType& vector()
  {
    CHECK_CRTP(this->as_imp().vector());
    return this->as_imp().vector();
  }

  const VectorType& vector() const
  {
    CHECK_CRTP(this->as_imp().vector());
    return this->as_imp().vector();
  }

  template< class S >
  ScalarType apply(const Stuff::LA::VectorInterface< S, ScalarType >& source) const
  {
    typedef typename S::derived_type SourceType;
    assemble();
    return vector().dot(static_cast< SourceType& >(source));
  }

  template< class S >
  ScalarType apply(const ConstDiscreteFunction< SpaceType, S >& source) const
  {
    assemble();
    assert(source.vector().size() == vector().size());
    return apply(source.vector());
  }
}; // class AssemblableFunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_INTERFACES_HH
