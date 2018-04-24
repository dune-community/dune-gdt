// This file is part of the dune-gdt project:
//   https://github.com/dune-community/dune-gdt
// Copyright 2010-2018 dune-gdt developers and contributors. All rights reserved.
// License: Dual licensed as BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
//      or  GPL-2.0+ (http://opensource.org/licenses/gpl-license)
//          with "runtime exception" (http://www.dune-project.org/license.html)
// Authors:
//   Felix Schindler (2015 - 2017)
//   Rene Milk       (2016, 2018)

#ifndef DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
#define DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH

#include <dune/xt/functions/reinterpret.hh>

#include "default.hh"

namespace Dune {
namespace GDT {


template <class DiscreteFunctionType>
class ReinterpretDiscreteFunction
    : public XT::Functions::ReinterpretFunction<DiscreteFunctionType,
                                                typename DiscreteFunctionType::SpaceType::GridLayerType>
{
  static_assert(is_const_discrete_function<DiscreteFunctionType>::value, "");
  typedef XT::Functions::ReinterpretFunction<DiscreteFunctionType,
                                             typename DiscreteFunctionType::SpaceType::GridLayerType>
      BaseType;

public:
  ReinterpretDiscreteFunction(const DiscreteFunctionType& source)
    : BaseType(source, source.space().grid_layer())
  {
  }
}; // class ReinterpretDiscreteFunction


} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_DISCRETEFUNCTION_REINTERPRET_HH
