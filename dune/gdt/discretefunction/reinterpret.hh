// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2016 dune-gdt developers and contributors. All rights reserved.
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
// Authors:
//   Felix Schindler (2015 - 2016)

#ifndef DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
#define DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH

#include <dune/stuff/functions/reinterpret.hh>

#include "default.hh"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType>
class ReinterpretDiscreteFunction
    : public Stuff::Functions::Reinterpret<DiscreteFunctionType, typename DiscreteFunctionType::SpaceType::GridViewType>
{
  static_assert(is_const_discrete_function<DiscreteFunctionType>::value, "");
  typedef Stuff::Functions::Reinterpret<DiscreteFunctionType, typename DiscreteFunctionType::SpaceType::GridViewType>
      BaseType;

public:
  ReinterpretDiscreteFunction(const DiscreteFunctionType& source)
    : BaseType(source, source.space().grid_view())
  {
  }
}; // class ReinterpretDiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
