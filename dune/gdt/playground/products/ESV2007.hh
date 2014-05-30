// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCTS_ESV2007_HH
#define DUNE_GDT_PRODUCTS_ESV2007_HH

#include <dune/gdt/playground/localevaluation/ESV2007.hh>
#include <dune/gdt/localoperator/codim0.hh>

#include "../../products/interfaces.hh"
#include "../../products/base.hh"

namespace Dune {
namespace GDT {
namespace Products {
namespace ESV2007 {


// forward, to be used in the traits
template< class GridViewImp,
          class DiffusionFactorType,
          class DiffusiveFluxType,
          class RangeImp,
          class SourceImp,
          class FieldImp = double,
          class DiffusionTensorType = void >
class DiffusiveFluxEstimate;


namespace internal {


template< class GridViewImp,
          class DiffusionFactorType,
          class DiffusiveFluxType,
          class RangeImp,
          class SourceImp,
          class FieldImp,
          class DiffusionTensorType = void >
class DiffusiveFluxEstimateTraits
{
public:
  typedef DiffusiveFluxEstimate
      < GridViewImp, DiffusionFactorType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp, DiffusionTensorType >
    derived_type;
  typedef GridViewImp GridViewType;
  typedef RangeImp    RangeType;
  typedef SourceImp   SourceType;
  typedef FieldImp    FieldType;
private:
  typedef LocalEvaluation::ESV2007::DiffusiveFluxEstimate
    < DiffusionFactorType, DiffusiveFluxType, DiffusionTensorType > LocalEvalautionType;
public:
  typedef LocalOperator::Codim0Integral< LocalEvalautionType > LocalOperatorType;
}; // class DiffusiveFluxEstimateTraits


template< class GridViewImp,
          class DiffusionType,
          class DiffusiveFluxType,
          class RangeImp,
          class SourceImp,
          class FieldImp >
class DiffusiveFluxEstimateTraits< GridViewImp, DiffusionType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp, void >
{
public:
  typedef DiffusiveFluxEstimate< GridViewImp, DiffusionType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp >
    derived_type;
  typedef GridViewImp GridViewType;
  typedef RangeImp    RangeType;
  typedef SourceImp   SourceType;
  typedef FieldImp    FieldType;
private:
  typedef LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionType, DiffusiveFluxType > LocalEvalautionType;
public:
  typedef LocalOperator::Codim0Integral< LocalEvalautionType > LocalOperatorType;
}; // class DiffusiveFluxEstimateTraits< ..., void >


} // namespace internal


template< class GridViewImp,
          class DiffusionType,
          class DiffusiveFluxType,
          class RangeImp,
          class SourceImp,
          class FieldImp >
class DiffusiveFluxEstimate< GridViewImp, DiffusionType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp, void >
  : public Products::LocalizableBase< internal::DiffusiveFluxEstimateTraits< GridViewImp,
                                                                             DiffusionType,
                                                                             DiffusiveFluxType,
                                                                             RangeImp,
                                                                             SourceImp,
                                                                             FieldImp > >
{
  typedef Products::LocalizableBase< internal::DiffusiveFluxEstimateTraits< GridViewImp,
                                                                            DiffusionType,
                                                                            DiffusiveFluxType,
                                                                            RangeImp,
                                                                            SourceImp,
                                                                            FieldImp > >
    BaseType;
public:
  typedef internal::DiffusiveFluxEstimateTraits
    < GridViewImp, DiffusionType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp > Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeType;
  using typename BaseType::SourceType;
  typedef typename Traits::LocalOperatorType  LocalOperatorType;

  DiffusiveFluxEstimate(const GridViewType& grd_vw,
                        const RangeType& rng,
                        const SourceType& src,
                        const DiffusionType& diffusion,
                        const DiffusiveFluxType& diffusive_flux,
                        const size_t over_integrate = 0)
    : BaseType(grd_vw, rng, src)
    , diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_, diffusive_flux_)
  {}

  ~DiffusiveFluxEstimate() {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return local_operator_;
  }

  const DiffusionType& diffusion_;
  const DiffusiveFluxType& diffusive_flux_;
  const LocalOperatorType local_operator_;
}; // class DiffusiveFluxEstimate< ..., void >


template< class GridViewImp,
          class DiffusionFactorType,
          class DiffusiveFluxType,
          class RangeImp,
          class SourceImp,
          class FieldImp,
          class DiffusionTensorType >
class DiffusiveFluxEstimate
  : public Products::LocalizableBase< internal::DiffusiveFluxEstimateTraits< GridViewImp,
                                                                             DiffusionFactorType,
                                                                             DiffusiveFluxType,
                                                                             RangeImp,
                                                                             SourceImp,
                                                                             FieldImp,
                                                                             DiffusionTensorType > >
{
  typedef Products::LocalizableBase< internal::DiffusiveFluxEstimateTraits< GridViewImp,
                                                                            DiffusionFactorType,
                                                                            DiffusiveFluxType,
                                                                            RangeImp,
                                                                            SourceImp,
                                                                            FieldImp,
                                                                            DiffusionTensorType > >
    BaseType;
public:
  typedef internal::DiffusiveFluxEstimateTraits
    < GridViewImp, DiffusionFactorType, DiffusiveFluxType, RangeImp, SourceImp, FieldImp, DiffusionTensorType > Traits;
  using typename BaseType::GridViewType;
  using typename BaseType::RangeType;
  using typename BaseType::SourceType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;

  DiffusiveFluxEstimate(const GridViewType& grd_vw,
                        const RangeType& rng,
                        const SourceType& src,
                        const DiffusionFactorType& diffusion_factor,
                        const DiffusionTensorType& diffusion_tensor,
                        const DiffusiveFluxType& diffusive_flux,
                        const size_t over_integrate = 0)
    : BaseType(grd_vw, rng, src)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_factor_, diffusion_tensor_, diffusive_flux_)
  {}

  ~DiffusiveFluxEstimate() {}

private:
  virtual const LocalOperatorType& local_operator() const DS_OVERRIDE DS_FINAL
  {
    return local_operator_;
  }

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const DiffusiveFluxType& diffusive_flux_;
  const LocalOperatorType local_operator_;
}; // class DiffusiveFluxEstimate


} // namespace ESV2007
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ESV2007_HH
