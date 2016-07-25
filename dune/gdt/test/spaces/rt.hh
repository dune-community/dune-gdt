// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2014 - 2016)

#ifndef DUNE_GDT_TEST_SPACES_RT_HH
#define DUNE_GDT_TEST_SPACES_RT_HH

#include <dune/common/unused.hh>

#include <dune/stuff/common/ranges.hh>

#include <dune/gdt/spaces/rt/interface.hh>

#include "base.hh"


template <class SpaceType>
class RT_Space : public LeafSpaceBase<SpaceType>
{
  template <class T, size_t d, size_t r, size_t rC>
  void matches_signature(const Dune::GDT::RtSpaceInterface<T, d, r, rC>& /*space*/)
  {
    static_assert(Dune::GDT::is_rt_space<SpaceType>::value, "");
    static_assert(std::is_same<typename SpaceType::Traits, T>::value, "");
    static_assert(d == SpaceType::dimDomain, "");
    static_assert(r == SpaceType::dimRange, "");
    static_assert(rC == SpaceType::dimRangeCols, "");
  }

public:
  virtual ~RT_Space()
  {
  }

  void matches_raviart_thomas_signature()
  {
    matches_signature(this->space_);
  }
};


template <class SpaceType>
class RT_2d_simplicial_Space : public RT_Space<SpaceType>
{
public:
  virtual ~RT_2d_simplicial_Space()
  {
  }

  void fulfills_raviart_thomas_2d_simplicial_interface()
  {
    for (const auto& entity : DSC::entityRange(this->space_.grid_view()))
      auto DUNE_UNUSED(local_DoF_indices) = this->space_.local_DoF_indices(entity);
  }
}; // class RT_2d_simplicial_Space


#endif // DUNE_GDT_TEST_SPACES_RT_HH
