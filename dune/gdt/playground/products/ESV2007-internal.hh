// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH

#include <dune/gdt/playground/localevaluation/ESV2007.hh>
#include <dune/gdt/localoperator/codim0.hh>

#include "../../products/base-internal.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace ESV2007 {
namespace internal {


template <class DiffusionFactorType, class DiffusiveFluxType, class GV, class FieldImp,
          class DiffusionTensorType = void>
class DiffusiveFluxEstimateBase : public LocalOperatorProviderBase<GV>
{
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, DiffusionFactorType>::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, DiffusionTensorType>::value,
                "DiffusionTensorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, DiffusiveFluxType>::value,
                "DiffusiveFluxType has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef DiffusiveFluxEstimateBase<DiffusionFactorType, DiffusiveFluxType, GV, FieldImp, DiffusionTensorType> ThisType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::
      Codim0Integral<LocalEvaluation::ESV2007::DiffusiveFluxEstimate<DiffusionFactorType, DiffusiveFluxType,
                                                                     DiffusionTensorType>> VolumeOperatorType;

  static const bool has_volume_operator = true;

  DiffusiveFluxEstimateBase(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                            const DiffusiveFluxType& diffusive_flux, const size_t over_integrate = 0)
    : volume_operator_(over_integrate, diffusion_factor, diffusion_tensor, diffusive_flux)
  {
  }

  DiffusiveFluxEstimateBase(const ThisType& other) = default;

  const VolumeOperatorType volume_operator_;
}; // class DiffusiveFluxEstimateBase


template <class DiffusionImp, class DiffusiveFluxImp, class GV, class FieldImp>
class DiffusiveFluxEstimateBase<DiffusionImp, DiffusiveFluxImp, GV, FieldImp, void>
    : public LocalOperatorProviderBase<GV>
{
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, DiffusionImp>::value,
                "DiffusionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::Tags::LocalizableFunction, DiffusiveFluxImp>::value,
                "DiffusiveFluxImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef DiffusiveFluxEstimateBase<DiffusionImp, DiffusiveFluxImp, GV, FieldImp, void> ThisType;

protected:
  typedef DiffusionImp DiffusionType;
  typedef DiffusiveFluxImp DiffusiveFluxType;

public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::ESV2007::DiffusiveFluxEstimate<DiffusionType,
                                                                                        DiffusiveFluxType>>
      VolumeOperatorType;

  static const bool has_volume_operator = true;

  DiffusiveFluxEstimateBase(const DiffusionType& diffusion, const DiffusiveFluxType& diffusive_flux,
                            const size_t over_integrate = 0)
    : volume_operator_(over_integrate, diffusion, diffusive_flux)
  {
  }

  DiffusiveFluxEstimateBase(const ThisType& other) = default;

  const VolumeOperatorType volume_operator_;
}; // class DiffusiveFluxEstimateBase


} // namespace internal
} // namespace ESV2007
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH
