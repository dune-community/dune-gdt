// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2017 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2014 - 2017)
//   Rene Milk       (2016)
//   Tobias Leibner  (2016)

#ifndef DUNE_GDT_TEST_SPACES_RT_HH
#define DUNE_GDT_TEST_SPACES_RT_HH

#include <dune/common/unused.hh>

#include <dune/xt/common/ranges.hh>

#include <dune/gdt/spaces/rt/interface.hh>

#include "base.hh"


template <class SpaceType>
class RT_Space : public SpaceBase<SpaceType>
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
    for (const auto& entity : elements(this->space_.grid_view()))
      auto local_DoF_indices DUNE_UNUSED = this->space_.local_DoF_indices(entity);
  }
}; // class RT_2d_simplicial_Space


#endif // DUNE_GDT_TEST_SPACES_RT_HH
