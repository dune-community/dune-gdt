// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ESV2007_HH
#define DUNE_GDT_PRODUCTS_ESV2007_HH

#include "ESV2007-internal.hh"
#include "../../products/base.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace ESV2007 {


// forward, to be used in the traits
template <class GridView, class DiffusionFactor, class DiffusiveFlux, class Range, class Source,
          class FieldType = double, class DiffusionTensor = void>
class DiffusiveFluxEstimate
    : public LocalizableBase<internal::DiffusiveFluxEstimateBase<DiffusionFactor, DiffusiveFlux, GridView, FieldType,
                                                                 DiffusionTensor>,
                             Range, Source>
{
  typedef LocalizableBase<internal::DiffusiveFluxEstimateBase<DiffusionFactor, DiffusiveFlux, GridView, FieldType,
                                                              DiffusionTensor>,
                          Range, Source> BaseType;

public:
  template <class... Args>
  DiffusiveFluxEstimate(Args&&... args)
    : BaseType(std::forward<Args>(args)...)
  {
  }
};


} // namespace ESV2007
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ESV2007_HH
