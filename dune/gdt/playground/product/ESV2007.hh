// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Albrecht
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PRODUCT_ESV2007_HH
#define DUNE_GDT_PRODUCT_ESV2007_HH

#include <type_traits>
#include <cmath>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/quadraturerules.hh>

#include <dune/stuff/functions/interfaces.hh>

#include <dune/gdt/playground/localevaluation/ESV2007.hh>
#include <dune/gdt/localoperator/codim0.hh>

#include "../../product/interfaces.hh"
#include "../../product/base.hh"

namespace Dune {
namespace GDT {
namespace Product {
namespace ESV2007 {


// forward, to be used in the traits
template <class GridViewImp, class FunctionImp>
class WeightedL2;


template <class GridViewImp, class FunctionImp>
class WeightedL2Traits
{
  static_assert(std::is_base_of<Stuff::IsLocalizableFunction, FunctionImp>::value,
                "FunctionImp has to be derived from Stuff::IsLocalizableFunction!");
  static_assert(std::is_same<typename GridViewImp::ctype, typename FunctionImp::DomainFieldType>::value,
                "Types do not match!");
  static_assert(GridViewImp::dimension == FunctionImp::dimDomain, "Dimensions do not match!");

public:
  typedef WeightedL2<GridViewImp, FunctionImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef FunctionImp FunctionType;
  typedef typename FunctionType::RangeFieldType FieldType;
}; // class WeightedL2Traits


template <class GridViewImp, class FunctionImp>
class WeightedL2 : public ProductInterface<WeightedL2Traits<GridViewImp, FunctionImp>>
{
public:
  typedef WeightedL2Traits<GridViewImp, FunctionImp> Traits;
  typedef typename Traits::GridViewType GridViewType;
  typedef typename Traits::FunctionType FunctionType;
  typedef typename Traits::FieldType FieldType;

  typedef typename GridViewType::template Codim<0>::Entity EntityType;
  typedef typename GridViewType::ctype DomainFieldType;
  static const unsigned int dimDomain = GridViewType::dimension;

  WeightedL2(const GridViewType& grid_view, const FunctionType& function,
             const FieldType poincare_constant = 1.0 / (M_PIl * M_PIl))
    : grid_view_(grid_view)
    , function_(function)
    , poincare_constant_(poincare_constant)
  {
    assert(poincare_constant_ > 0);
  }

  template <class RR, int rRR, int rCR, class RS, int rRS, int rCS>
  FieldType apply2(
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RR, rRR, rCR>& /*range*/,
      const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, RS, rRS, rCS>& /*source*/) const
  {
    static_assert((Dune::AlwaysFalse<RR>::value), "Not implemented for this combination!");
  }

  FieldType
  apply2(const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1, 1>& range,
         const Stuff::LocalizableFunctionInterface<EntityType, DomainFieldType, dimDomain, FieldType, 1, 1>& source,
         const size_t over_integrate = 0) const
  {
    typedef typename FunctionType::RangeType FunctionRangeType;
    FunctionRangeType function_value(0);
    FieldType min_eigen_value(std::numeric_limits<FieldType>::max());
    FieldVector<FieldType, 1> range_value(0);
    FieldVector<FieldType, 1> source_value(0);
    FieldType ret(0);
    // walk the grid
    const auto entity_it_end = grid_view_.template end<0>();
    for (auto entity_it = grid_view_.template begin<0>(); entity_it != entity_it_end; ++entity_it) {
      const auto& entity = *entity_it;
      // get local functions
      const auto local_function = function_.local_function(entity);
      const auto local_range    = range.local_function(entity);
      const auto local_source   = source.local_function(entity);
      // do a volume quadrature
      FieldType integral(0);
      const size_t integrand_order =
          std::max(local_function->order(), std::max(local_range->order(), local_source->order())) + over_integrate;
      assert(integrand_order < std::numeric_limits<int>::max());
      const auto& quadrature       = QuadratureRules<DomainFieldType, dimDomain>::rule(entity.type(), int(integrand_order));
      const auto quadrature_it_end = quadrature.end();
      for (auto quadrature_it = quadrature.begin(); quadrature_it != quadrature_it_end; ++quadrature_it) {
        const FieldVector<DomainFieldType, dimDomain> xx = quadrature_it->position();
        const double integration_element                 = entity.geometry().integrationElement(xx);
        const double quadrature_weight                   = quadrature_it->weight();
        // find minimum of weighting function
        local_function->evaluate(xx, function_value);
        min_eigen_value = std::min(min_eigen_value, compute_min_eigen_value_of_(function_value));
        // compute integral
        local_range->evaluate(xx, range_value);
        local_source->evaluate(xx, source_value);
        integral += integration_element * quadrature_weight * (range_value * source_value);
      } // do a volume quadrature
      assert(min_eigen_value > 0);
      const DomainFieldType hh         = compute_diameter_of_(entity);
      const FieldType weighting_factor = (poincare_constant_ * hh * hh) / min_eigen_value;
      ret += weighting_factor * integral;
    } // walk the grid
    return ret;
  } // ... apply2(...)

private:
  DomainFieldType compute_diameter_of_(const EntityType& entity) const
  {
    DomainFieldType ret(0);
    for (int cc = 0; cc < entity.template count<dimDomain>(); ++cc) {
      const auto vertex = entity.template subEntity<dimDomain>(cc)->geometry().center();
      for (int dd = cc + 1; dd < entity.template count<dimDomain>(); ++dd) {
        const auto other_vertex = entity.template subEntity<dimDomain>(dd)->geometry().center();
        const auto diff         = vertex - other_vertex;
        ret                     = std::max(ret, diff.two_norm());
      }
    }
    return ret;
  } // ... compute_diameter(...)

  FieldType compute_min_eigen_value_of_(const FieldVector<FieldType, 1>& value) const
  {
    return value[0];
  }

  template <class F>
  FieldType compute_min_eigen_value(const FieldMatrix<F, dimDomain, dimDomain>& /*value*/) const
  {
    // the template is only here for the static assert, replace by FieldType upon implementing
    static_assert(AlwaysFalse<F>::value, "Not yet implemented!");
  }

  const GridViewType& grid_view_;
  const FunctionType& function_;
  const FieldType poincare_constant_;
}; // class WeightedL2


// forward, to be used in the traits
template <class GridViewImp, class DiffusionImp, class DiffusiveFluxImp, class RangeImp, class SourceImp,
          class FieldImp = double>
class DiffusiveFluxEstimate;


template <class GridViewImp, class DiffusionImp, class DiffusiveFluxImp, class RangeImp, class SourceImp,
          class FieldImp>
class DiffusiveFluxEstimateTraits
{
public:
  typedef DiffusiveFluxEstimate<GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp>
      derived_type;
  typedef GridViewImp GridViewType;
  typedef DiffusionImp DiffusionType;
  typedef DiffusiveFluxImp DiffusiveFluxType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;

private:
  typedef LocalEvaluation::ESV2007::DiffusiveFluxEstimate<DiffusionType, DiffusiveFluxType> LocalEvalautionType;

public:
  typedef LocalOperator::Codim0Integral<LocalEvalautionType> LocalOperatorType;
}; // class DiffusiveFluxEstimateTraits


template <class GridViewImp, class DiffusionImp, class DiffusiveFluxImp, class RangeImp, class SourceImp,
          class FieldImp>
class DiffusiveFluxEstimate
    : public Product::LocalizableBase<DiffusiveFluxEstimateTraits<GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp,
                                                                  SourceImp, FieldImp>>
{
  typedef Product::LocalizableBase<DiffusiveFluxEstimateTraits<GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp,
                                                               SourceImp, FieldImp>> BaseType;

public:
  typedef DiffusiveFluxEstimateTraits<GridViewImp, DiffusionImp, DiffusiveFluxImp, RangeImp, SourceImp, FieldImp>
      Traits;
  using typename BaseType::GridViewType;
  typedef typename Traits::DiffusionType DiffusionType;
  typedef typename Traits::DiffusiveFluxType DiffusiveFluxType;
  using typename BaseType::RangeType;
  using typename BaseType::SourceType;
  typedef typename Traits::LocalOperatorType LocalOperatorType;


  DiffusiveFluxEstimate(const GridViewType& grid_view, const RangeType& range, const SourceType& source,
                        const DiffusionType& diffusion, const DiffusiveFluxType& diffusive_flux,
                        const size_t over_integrate = 0)
    : BaseType(grid_view, range, source)
    , diffusion_(diffusion)
    , diffusive_flux_(diffusive_flux)
    , local_operator_(over_integrate, diffusion_, diffusive_flux_)
  {
  }

  ~DiffusiveFluxEstimate()
  {
  }

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
} // namespace Product
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PRODUCT_ESV2007_HH
