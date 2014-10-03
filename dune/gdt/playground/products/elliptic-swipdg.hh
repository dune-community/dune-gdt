// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_PRODUCTS_BASE_HH
#define DUNE_GDT_PLAYGROUND_PRODUCTS_BASE_HH

#include <type_traits>

#include <dune/grid/common/gridview.hh>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/common/tmp-storage.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/grid/walker.hh>
#include <dune/stuff/la/container/pattern.hh>

#include <dune/gdt/assembler/system.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/playground/localevaluation/swipdg.hh>
#include <dune/gdt/products/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Products {


// forward, needed for the traits
template <class GridViewImp, class DiffusionFactorImp, class RangeImp, class SourceImp, class FieldImp,
          class DiffusionTensorImp>
class EllipticSWIPDGPenaltyLocalizable;


namespace internal {


template <class GridViewImp, class DiffusionFactorImp, class RangeImp, class SourceImp, class FieldImp,
          class DiffusionTensorImp>
class EllipticSWIPDGPenaltyLocalizableTraits
{
public:
  typedef EllipticSWIPDGPenaltyLocalizable<GridViewImp, DiffusionFactorImp, RangeImp, SourceImp, FieldImp,
                                           DiffusionTensorImp> derived_type;
  typedef GridViewImp GridViewType;
  typedef RangeImp RangeType;
  typedef SourceImp SourceType;
  typedef FieldImp FieldType;
}; // class EllipticSWIPDGPenaltyLocalizableTraits


} // namespace internal


template <class GridViewImp, class DiffusionFactorImp, class RangeImp, class SourceImp = RangeImp,
          class FieldImp = double, class DiffusionTensorImp = void>
class EllipticSWIPDGPenaltyLocalizable
    : public LocalizableProductInterface<internal::EllipticSWIPDGPenaltyLocalizableTraits<GridViewImp,
                                                                                          DiffusionFactorImp, RangeImp,
                                                                                          SourceImp, FieldImp,
                                                                                          DiffusionTensorImp>>,
      public Stuff::Grid::Functor::Codim1<GridViewImp>
{
  typedef LocalizableProductInterface<internal::EllipticSWIPDGPenaltyLocalizableTraits<GridViewImp, DiffusionFactorImp,
                                                                                       RangeImp, SourceImp, FieldImp,
                                                                                       DiffusionTensorImp>>
      ProductBaseType;
  typedef Stuff::Grid::Functor::Codim1<GridViewImp> FunctorBaseType;

public:
  typedef internal::EllipticSWIPDGPenaltyLocalizableTraits<GridViewImp, DiffusionFactorImp, RangeImp, SourceImp,
                                                           FieldImp, DiffusionTensorImp> Traits;

  typedef typename FunctorBaseType::GridViewType GridViewType;
  typedef typename ProductBaseType::RangeType RangeType;
  typedef typename ProductBaseType::SourceType SourceType;
  typedef typename ProductBaseType::FieldType FieldType;
  typedef typename FunctorBaseType::EntityType EntityType;
  typedef typename FunctorBaseType::IntersectionType IntersectionType;
  typedef DiffusionFactorImp DiffusionFactorType;
  typedef DiffusionTensorImp DiffusionTensorType;

  typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDG::InnerPenalty<DiffusionFactorType,
                                                                                      DiffusionTensorType>>
      CouplingOperatorType;
  typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SWIPDG::BoundaryLHSPenalty<DiffusionFactorType,
                                                                                            DiffusionTensorType>>
      BoundaryOperatorType;

private:
  typedef DSC::TmpMatricesStorage<FieldType> TmpMatricesProviderType;

public:
  EllipticSWIPDGPenaltyLocalizable(const GridViewType& grd_vw, const RangeType& rng, const SourceType& src,
                                   const DiffusionFactorType& diffusion_factor,
                                   const DiffusionTensorType& diffusion_tensor, const size_t over_integrate = 0)
    : grid_view_(grd_vw)
    , range_(rng)
    , source_(src)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , coupling_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
    , boundary_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
    , inner_intersections_()
    , boundary_intersections_()
    , tmp_storage_(nullptr)
    , prepared_(false)
    , finalized_(false)
    , result_(0)
    , finalized_result_(0)
  {
  }

  virtual ~EllipticSWIPDGPenaltyLocalizable()
  {
  }

  const GridViewType& grid_view() const
  {
    return grid_view_;
  }

  const RangeType& range() const
  {
    return range_;
  }

  const SourceType& source() const
  {
    return source_;
  }

  virtual void prepare() DS_OVERRIDE
  {
    if (!prepared_) {
      tmp_storage_ = std::unique_ptr<TmpMatricesProviderType>(new TmpMatricesProviderType(
          {4, std::max(coupling_operator_.numTmpObjectsRequired(), boundary_operator_.numTmpObjectsRequired())}, 1, 1));
      result_   = FieldType(0.0);
      prepared_ = true;
    }
  } // ... prepare()

  FieldType compute_locally(const IntersectionType& intersection, const EntityType& inside_entity,
                            const EntityType& outside_entity) const
  {
    assert(prepared_);
    assert(tmp_storage_);
    auto& tmp_storage = tmp_storage_->matrices();
    assert(tmp_storage.size() >= 2);
    assert(tmp_storage[0].size() >= 4);
    auto& local_operator_result_en_en = tmp_storage[0][0];
    auto& local_operator_result_ne_ne = tmp_storage[0][1];
    auto& local_operator_result_en_ne = tmp_storage[0][2];
    auto& local_operator_result_ne_en = tmp_storage[0][3];
    auto& tmp_matrices                = tmp_storage[1];
    // get the local functions
    const auto local_source_ptr_en = this->source().local_function(inside_entity);
    const auto local_source_ptr_ne = this->source().local_function(outside_entity);
    const auto local_range_ptr_en  = this->range().local_function(inside_entity);
    const auto local_range_ptr_ne  = this->range().local_function(outside_entity);
    // apply local operator
    if (inner_intersections_.apply_on(grid_view_, intersection)) {
      coupling_operator_.apply(*local_range_ptr_en,
                               *local_source_ptr_en,
                               *local_source_ptr_ne,
                               *local_range_ptr_ne,
                               intersection,
                               local_operator_result_en_en,
                               local_operator_result_ne_ne,
                               local_operator_result_en_ne,
                               local_operator_result_ne_en,
                               tmp_matrices);
      assert(local_operator_result_en_en.rows() == 1);
      assert(local_operator_result_en_en.cols() == 1);
      assert(local_operator_result_ne_ne.rows() == 1);
      assert(local_operator_result_ne_ne.cols() == 1);
      assert(local_operator_result_en_ne.rows() == 1);
      assert(local_operator_result_en_ne.cols() == 1);
      assert(local_operator_result_ne_en.rows() == 1);
      assert(local_operator_result_ne_en.cols() == 1);
      return local_operator_result_en_en[0][0] + local_operator_result_ne_ne[0][0] + local_operator_result_en_ne[0][0]
             + local_operator_result_ne_en[0][0];
    }
    if (boundary_intersections_.apply_on(grid_view_, intersection)) {
      boundary_operator_.apply(
          *local_range_ptr_en, *local_source_ptr_en, intersection, local_operator_result_en_en, tmp_matrices);
      assert(local_operator_result_en_en.rows() == 1);
      assert(local_operator_result_en_en.cols() == 1);
      return local_operator_result_en_en[0][0];
    }
  } // ... compute_locally(...)

  virtual void apply_local(const IntersectionType& intersection, const EntityType& inside_entity,
                           const EntityType& outside_entity) DS_OVERRIDE
  {
    *result_ += compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() DS_OVERRIDE
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_        = true;
    }
  } // ... finalize(...)

  FieldType apply2()
  {
    if (!finalized_) {
      Stuff::Grid::Walker<GridViewType> grid_walker(grid_view_);
      grid_walker.add(*this);
      grid_walker.walk();
    }
    return finalized_result_;
  } // ... apply2(...)

private:
  const GridViewType& grid_view_;
  const RangeType& range_;
  const SourceType& source_;
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const CouplingOperatorType coupling_operator_;
  const BoundaryOperatorType boundary_operator_;
  const DSG::ApplyOn::InnerIntersectionsPrimally<GridViewType> inner_intersections_;
  const DSG::ApplyOn::BoundaryIntersections<GridViewType> boundary_intersections_;
  std::unique_ptr<TmpMatricesProviderType> tmp_storage_;
  bool prepared_;
  bool finalized_;
  DS::PerThreadValue<FieldType> result_;
  FieldType finalized_result_;
}; // class EllipticSWIPDGPenaltyLocalizable


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_PRODUCTS_BASE_HH
