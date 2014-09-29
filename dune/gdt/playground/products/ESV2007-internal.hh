// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH
#define DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH

#include <dune/gdt/playground/localevaluation/ESV2007.hh>
#include <dune/gdt/localoperator/codim0.hh>

#include "../../products/base.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace ESV2007 {
namespace internal {


template< class DiffusionFactorType, class DiffusiveFluxType, class GV, class FieldImp, class DiffusionTensorType = void >
class DiffusiveFluxEstimateBase
{
  static_assert(std::is_base_of< Stuff::Tags::LocalizableFunction, DiffusionFactorType >::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of< Stuff::Tags::LocalizableFunction, DiffusionTensorType >::value,
                "DiffusionTensorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of< Stuff::Tags::LocalizableFunction, DiffusiveFluxType >::value,
                "DiffusiveFluxType has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef DiffusiveFluxEstimateBase< DiffusionFactorType, DiffusiveFluxType, GV, FieldImp, DiffusionTensorType >
      ThisType;
public:
  typedef GV       GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral
      < LocalEvaluation::ESV2007::DiffusiveFluxEstimate
          < DiffusionFactorType, DiffusiveFluxType, DiffusionTensorType > > LocalOperatorType;

  DiffusiveFluxEstimateBase(const DiffusionFactorType& diffusion_factor,
                            const DiffusionTensorType& diffusion_tensor,
                            const DiffusiveFluxType& diffusive_flux,
                            const size_t over_integrate = 0)
    : diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_factor_, diffusion_tensor_, diffusive_flux_)
  {}

  DiffusiveFluxEstimateBase(const ThisType& other) = default;

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
protected:
  const LocalOperatorType local_operator_;
}; // class DiffusiveFluxEstimateBase


template< class DiffusionImp, class DiffusiveFluxImp, class GV, class FieldImp >
class DiffusiveFluxEstimateBase< DiffusionImp, DiffusiveFluxImp, GV, FieldImp, void >
{
  static_assert(std::is_base_of< Stuff::Tags::LocalizableFunction, DiffusionImp >::value,
                "DiffusionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of< Stuff::Tags::LocalizableFunction, DiffusiveFluxImp >::value,
                "DiffusiveFluxImp has to be derived from Stuff::LocalizableFunctionInterface!");
  typedef DiffusiveFluxEstimateBase< DiffusionImp, DiffusiveFluxImp, GV, FieldImp, void > ThisType;
protected:
  typedef DiffusionImp     DiffusionType;
  typedef DiffusiveFluxImp DiffusiveFluxType;
public:
  typedef GV GridViewType;
  typedef FieldImp FieldType;
  typedef LocalOperator::Codim0Integral
      < LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionType, DiffusiveFluxType > > LocalOperatorType;

  DiffusiveFluxEstimateBase(const DiffusionType& diffusion,
                            const DiffusiveFluxType& diffusive_flux,
                            const size_t over_integrate = 0)
    : diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_, diffusive_flux_)
  {}

  DiffusiveFluxEstimateBase(const ThisType& other) = default;

private:
  const DiffusionType& diffusion_;
  const DiffusiveFluxType& diffusive_flux_;
protected:
  const LocalOperatorType local_operator_;
}; // class DiffusiveFluxEstimateBase


} // namespace internal
} // namespace ESV2007
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ESV2007_INTERNAL_HH
