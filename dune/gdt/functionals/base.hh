// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_FUNCTIONALS_BASE_HH
#define DUNE_GDT_FUNCTIONALS_BASE_HH

#include <dune/stuff/la/solver.hh>

#include "interfaces.hh"

namespace Dune {
namespace GDT {
namespace Functionals {


template <class Traits>
class VectorBased : public AssemblableFunctionalInterface<Traits>
{
  typedef AssemblableFunctionalInterface<Traits> BaseType;

public:
  using typename BaseType::GridViewType;
  using typename BaseType::SpaceType;
  using typename BaseType::VectorType;
  using typename BaseType::ScalarType;

  VectorBased(VectorType& vec, const SpaceType& spc, const GridViewType& grd_vw)
    : vector_(vec)
    , space_(spc)
    , grid_view_(grd_vw)
  {
  }

  VectorBased(VectorType& vec, const SpaceType& spc)
    : vector_(vec)
    , space_(spc)
    , grid_view_(spc.grid_view())
  {
  }

  virtual ~VectorBased()
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const SpaceType& space() const
  {
    return space_;
  }

  VectorType& vector()
  {
    return vector_;
  }

  const VectorType& vector() const
  {
    return vector_;
  }

  virtual void assemble() = 0;

  template <class S>
  ScalarType apply(const Stuff::LA::VectorInterface<S>& source) const
  {
    assemble();
    return vector_.dot(source.as_imp());
  } // ... apply(...)

private:
  VectorType& vector_;
  const SpaceType& space_;
  const GridViewType& grid_view_;
}; // class VectorBased


// forward, to be used in the traits
template <class ImpTraits>
class AssemblableFaceBase;


template <class ImpTraits>
class AssemblableFaceBaseTraits
{
public:
  typedef typename ImpTraits::derived_type derived_type;
  typedef typename ImpTraits::GridViewType GridViewType;
  typedef typename ImpTraits::SpaceType SpaceType;
  typedef typename ImpTraits::VectorType VectorType;
  typedef typename ImpTraits::ScalarType ScalarType;

private:
  static_assert(std::is_base_of<Stuff::LA::VectorInterface<typename VectorType::Traits>, VectorType>::value,
                "VectorType has to be derived from Stuff::LA::VectorInterface!");
}; // class AssemblableFaceBaseTraits


} // namespace Functionals
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_FUNCTIONALS_BASE_HH
