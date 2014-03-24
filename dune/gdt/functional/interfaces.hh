// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONAL_INTERFACES_HH
#define DUNE_GDT_FUNCTIONAL_INTERFACES_HH

#include <dune/stuff/common/crtp.hh>

namespace Dune {
namespace GDT {


template <class Traits>
class FunctionalInterface : public Stuff::CRTPInterface<FunctionalInterface<Traits>, Traits>
{
  typedef typename Traits::derived_type derived_type;

public:
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::ScalarType ScalarType;

  template <class SourceType>
  ScalarType apply(const SourceType& source) const
  {
    CHECK_CRTP(this->as_imp(*this).apply(source));
    return this->as_imp(*this).apply(source);
  } // apply(...)
}; // class FunctionalInterface


template <class Traits>
class AssemblableFunctionalInterface : protected Stuff::CRTPInterface<AssemblableFunctionalInterface<Traits>, Traits>
{
public:
  typedef typename Traits::derived_type derived_type;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::SpaceType SpaceType;
  typedef typename Traits::VectorType VectorType;
  typedef typename Traits::ScalarType ScalarType;

private:
  static_assert(std::is_base_of<SpaceInterface<typename SpaceType::Traits>, SpaceType>::value,
                "SpaceType has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<Stuff::LA::VectorInterface<typename VectorType::Traits>, VectorType>::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");

public:
  const GridViewType& grid_view() const
  {
    CHECK_CRTP(this->as_imp(*this).grid_view());
    return this->as_imp(*this).grid_view();
  }

  const SpaceType& space() const
  {
    CHECK_CRTP(this->as_imp(*this).space());
    return this->as_imp(*this).space();
  }

  void assemble()
  {
    CHECK_AND_CALL_CRTP(this->as_imp(*this).assemble());
  }

  VectorType& vector()
  {
    CHECK_CRTP(this->as_imp(*this).vector());
    return this->as_imp(*this).vector();
  }

  const VectorType& vector() const
  {
    CHECK_CRTP(this->as_imp(*this).vector());
    return this->as_imp(*this).vector();
  }

  template <class S>
  ScalarType apply(const Stuff::LA::VectorInterface<S>& source) const
  {
    CHECK_CRTP(this->as_imp(*this).apply(source));
    return this->as_imp(*this).apply(source);
  }

  template <class S>
  ScalarType apply(const ConstDiscreteFunction<SpaceType, S>& source) const
  {
    return apply(source.vector());
  }
}; // class AssemblableFunctionalInterface


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONAL_INTERFACES_HH
