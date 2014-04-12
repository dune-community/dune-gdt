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
template< class GridViewImp, class DiffusionImp, class DiffusiveFluxImp
        , class RangeImp, class SourceImp, class FieldImp = double >
class DiffusiveFluxEstimate;


template< class GridViewImp, class DiffusionImp, class DiffusiveFluxImp
        , class RangeImp, class SourceImp, class FieldImp >
class DiffusiveFluxEstimateTraits
{
public:
  typedef DiffusiveFluxEstimate< GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp >
    derived_type;
  typedef GridViewImp GridViewType;
  typedef DiffusionImp DiffusionType;
  typedef DiffusiveFluxImp DiffusiveFluxType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;
private:
  typedef LocalEvaluation::ESV2007::DiffusiveFluxEstimate< DiffusionType, DiffusiveFluxType > LocalEvalautionType;
public:
  typedef LocalOperator::Codim0Integral< LocalEvalautionType > LocalOperatorType;
}; // class DiffusiveFluxEstimateTraits


template< class GridViewImp, class DiffusionImp, class DiffusiveFluxImp
        , class RangeImp, class SourceImp, class FieldImp >
class DiffusiveFluxEstimate
  : public Products::LocalizableBase< DiffusiveFluxEstimateTraits< GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp > >
{
  typedef Products::LocalizableBase< DiffusiveFluxEstimateTraits< GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp > >
    BaseType;
public:
  typedef DiffusiveFluxEstimateTraits< GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp >
    Traits;
  using typename BaseType::GridViewType;
  typedef typename Traits::DiffusionType      DiffusionType;
  typedef typename Traits::DiffusiveFluxType  DiffusiveFluxType;
  using typename BaseType::RangeType;
  using typename BaseType::SourceType;
  typedef typename Traits::LocalOperatorType  LocalOperatorType;


  DiffusiveFluxEstimate(const GridViewType& grid_view, const RangeType& range, const SourceType& source,
                        const DiffusionType& diffusion, const DiffusiveFluxType& diffusive_flux,
                        const size_t over_integrate = 0)
    : BaseType(grid_view, range, source)
    , diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_, diffusive_flux_)
  {}

  ~DiffusiveFluxEstimate() {}

private:
  virtual const LocalOperatorType& local_operator() const
  {
    return local_operator_;
  }

  const DiffusionType& diffusion_;
  const DiffusiveFluxType& diffusive_flux_;
  const LocalOperatorType local_operator_;
}; // class DiffusiveFluxEstimate


} // namespace ESV2007
} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCTS_ESV2007_HH
