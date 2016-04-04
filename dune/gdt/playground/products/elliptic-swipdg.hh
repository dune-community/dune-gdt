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
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/products/interfaces.hh>
#include <dune/gdt/spaces/interface.hh>

namespace Dune {
namespace GDT {
namespace Products {


// forward, needed for the traits
template< class GridViewImp,
          class DiffusionFactorImp,
          class RangeImp,
          class SourceImp,
          class FieldImp,
          class DiffusionTensorImp >
class EllipticSWIPDGPenaltyLocalizable;


template< class DiffusionFactorImp,
          class MatrixImp,
          class RangeSpaceImp,
          class GridViewImp,
          class SourceSpaceImp,
          class DiffusionTensorImp >
class EllipticSWIPDGPenaltyAssemblable;


namespace internal {


template< class GridViewImp,
          class DiffusionFactorImp,
          class RangeImp,
          class SourceImp,
          class FieldImp,
          class DiffusionTensorImp >
class EllipticSWIPDGPenaltyLocalizableTraits
{
public:
  typedef EllipticSWIPDGPenaltyLocalizable
      < GridViewImp, DiffusionFactorImp, RangeImp, SourceImp, FieldImp, DiffusionTensorImp > derived_type;
  typedef GridViewImp GridViewType;
  typedef RangeImp    RangeType;
  typedef SourceImp   SourceType;
  typedef FieldImp    FieldType;
}; // class EllipticSWIPDGPenaltyLocalizableTraits


template< class DiffusionFactorImp,
          class MatrixImp,
          class RangeSpaceImp,
          class GridViewImp,
          class SourceSpaceImp,
          class DiffusionTensorImp >
class EllipticSWIPDGPenaltyAssemblableTraits
{
public:
  typedef EllipticSWIPDGPenaltyAssemblable
      < DiffusionFactorImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, DiffusionTensorImp > derived_type;
  typedef GridViewImp    GridViewType;
  typedef RangeSpaceImp  RangeSpaceType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef MatrixImp      MatrixType;
}; // class EllipticSWIPDGPenaltyAssemblableTraits


} // namespace internal


template< class GridViewImp,
          class DiffusionFactorImp,
          class RangeImp,
          class SourceImp = RangeImp,
          class FieldImp = double,
          class DiffusionTensorImp = void >
class EllipticSWIPDGPenaltyLocalizable
  : public LocalizableProductInterface< internal::EllipticSWIPDGPenaltyLocalizableTraits< GridViewImp, DiffusionFactorImp, RangeImp, SourceImp, FieldImp, DiffusionTensorImp > >
  , public Stuff::Grid::Functor::Codim1< GridViewImp >
{
  typedef LocalizableProductInterface< internal::EllipticSWIPDGPenaltyLocalizableTraits
      < GridViewImp, DiffusionFactorImp, RangeImp, SourceImp, FieldImp, DiffusionTensorImp > > ProductBaseType;
  typedef Stuff::Grid::Functor::Codim1< GridViewImp >                                          FunctorBaseType;
public:
  typedef internal::EllipticSWIPDGPenaltyLocalizableTraits
      < GridViewImp, DiffusionFactorImp, RangeImp, SourceImp, FieldImp, DiffusionTensorImp >   Traits;

  typedef typename FunctorBaseType::GridViewType     GridViewType;
  typedef typename ProductBaseType::RangeType        RangeType;
  typedef typename ProductBaseType::SourceType       SourceType;
  typedef typename ProductBaseType::FieldType        FieldType;
  typedef typename FunctorBaseType::EntityType       EntityType;
  typedef typename FunctorBaseType::IntersectionType IntersectionType;
  typedef DiffusionFactorImp DiffusionFactorType;
  typedef DiffusionTensorImp DiffusionTensorType;

  typedef LocalOperator::Codim1CouplingIntegral
      < LocalEvaluation::SWIPDG::InnerPenalty< DiffusionFactorType, DiffusionTensorType > >       CouplingOperatorType;
  typedef LocalOperator::Codim1BoundaryIntegral
      < LocalEvaluation::SWIPDG::BoundaryLHSPenalty< DiffusionFactorType, DiffusionTensorType > > BoundaryOperatorType;

private:
  typedef DSC::TmpMatricesStorage< FieldType > TmpMatricesProviderType;

public:
  EllipticSWIPDGPenaltyLocalizable(const GridViewType& grd_vw,
                        const RangeType& rng,
                        const SourceType& src,
                        const DiffusionFactorType& diffusion_factor,
                        const DiffusionTensorType& diffusion_tensor,
                        const size_t over_integrate = 0)
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
  {}

  virtual ~EllipticSWIPDGPenaltyLocalizable() {}

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

  virtual void prepare() override
  {
    if (!prepared_) {
      tmp_storage_ = std::unique_ptr< TmpMatricesProviderType >(new TmpMatricesProviderType({4,
                                                                  std::max(coupling_operator_.numTmpObjectsRequired(),
                                                                           boundary_operator_.numTmpObjectsRequired())},
                                                                 1, 1));
      result_ = FieldType(0.0);
      prepared_ = true;
    }
  } // ... prepare()

  FieldType compute_locally(const IntersectionType& intersection,
                            const EntityType& inside_entity,
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
    auto& tmp_matrices = tmp_storage[1];
    // get the local functions
    const auto local_source_ptr_en = this->source().local_function(inside_entity);
    const auto local_source_ptr_ne = this->source().local_function(outside_entity);
    const auto local_range_ptr_en = this->range().local_function(inside_entity);
    const auto local_range_ptr_ne = this->range().local_function(outside_entity);
    // apply local operator
    FieldType ret = 0;
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
    ret += local_operator_result_en_en[0][0]
         + local_operator_result_ne_ne[0][0]
         + local_operator_result_en_ne[0][0]
         + local_operator_result_ne_en[0][0];
    }
    if (boundary_intersections_.apply_on(grid_view_, intersection)) {
      boundary_operator_.apply(*local_range_ptr_en,
                               *local_source_ptr_en,
                               intersection,
                               local_operator_result_en_en,
                               tmp_matrices);
      assert(local_operator_result_en_en.rows() == 1);
      assert(local_operator_result_en_en.cols() == 1);
      ret += local_operator_result_en_en[0][0];
    }
    return ret;
  } // ... compute_locally(...)

  virtual void apply_local(const IntersectionType& intersection,
                           const EntityType& inside_entity,
                           const EntityType& outside_entity) override
  {
    *result_ += compute_locally(intersection, inside_entity, outside_entity);
  }

  virtual void finalize() override
  {
    if (!finalized_) {
      finalized_result_ = result_.sum();
      finalized_result_ = grid_view_.comm().sum(finalized_result_);
      finalized_ = true;
    }
  } // ... finalize(...)

  FieldType apply2()
  {
    if (!finalized_) {
      Stuff::Grid::Walker< GridViewType > grid_walker(grid_view_);
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
  const DSG::ApplyOn::InnerIntersectionsPrimally< GridViewType > inner_intersections_;
  const DSG::ApplyOn::BoundaryIntersections< GridViewType > boundary_intersections_;
  std::unique_ptr< TmpMatricesProviderType > tmp_storage_;
  bool prepared_;
  bool finalized_;
  DS::PerThreadValue< FieldType > result_;
  FieldType finalized_result_;
}; // class EllipticSWIPDGPenaltyLocalizable


template< class DiffusionFactorImp,
          class MatrixImp,
          class RangeSpaceImp,
          class GridViewImp = typename RangeSpaceImp::GridViewType,
          class SourceSpaceImp = RangeSpaceImp,
          class DiffusionTensorImp = void >
class EllipticSWIPDGPenaltyAssemblable
  : DSC::StorageProvider< MatrixImp >
  , public AssemblableProductInterface< internal::EllipticSWIPDGPenaltyAssemblableTraits
        < DiffusionFactorImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, DiffusionTensorImp > >
  , public SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >
{
  typedef DSC::StorageProvider< MatrixImp >                                                 StorageBaseType;
  typedef AssemblableProductInterface< internal::EllipticSWIPDGPenaltyAssemblableTraits
      < DiffusionFactorImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, DiffusionTensorImp > >
                                                                                            ProductBaseType;
  typedef SystemAssembler < RangeSpaceImp, GridViewImp, SourceSpaceImp >                    AssemblerBaseType;
public:
  typedef internal::EllipticSWIPDGPenaltyAssemblableTraits
      < DiffusionFactorImp, MatrixImp, RangeSpaceImp, GridViewImp, SourceSpaceImp, DiffusionTensorImp > Traits;
  typedef typename AssemblerBaseType::GridViewType     GridViewType;
  typedef typename ProductBaseType::RangeSpaceType     RangeSpaceType;
  typedef typename ProductBaseType::SourceSpaceType    SourceSpaceType;
  typedef typename ProductBaseType::FieldType          FieldType;
  typedef typename ProductBaseType::MatrixType         MatrixType;
  typedef typename AssemblerBaseType::EntityType       EntityType;
  typedef typename AssemblerBaseType::IntersectionType IntersectionType;
  typedef DiffusionFactorImp DiffusionFactorType;
  typedef DiffusionTensorImp DiffusionTensorType;

  typedef LocalOperator::Codim1CouplingIntegral
      < LocalEvaluation::SWIPDG::InnerPenalty< DiffusionFactorType, DiffusionTensorType > >       CouplingOperatorType;
  typedef LocalOperator::Codim1BoundaryIntegral
      < LocalEvaluation::SWIPDG::BoundaryLHSPenalty< DiffusionFactorType, DiffusionTensorType > > BoundaryOperatorType;
  typedef LocalAssembler::Codim1CouplingMatrix< CouplingOperatorType >                            CouplingAssemblerType;
  typedef LocalAssembler::Codim1BoundaryMatrix< BoundaryOperatorType >                            BoundaryAssemblerType;

  using ProductBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return range_space.compute_face_and_volume_pattern(grid_view, source_space);
  }

  EllipticSWIPDGPenaltyAssemblable(const RangeSpaceType& rng_spc,
                                   const GridViewType& grd_vw,
                                   const SourceSpaceType& src_spc,
                                   const DiffusionFactorImp& diffusion_factor,
                                   const DiffusionTensorImp& diffusion_tensor,
                                   const size_t over_integrate = 0)
    : StorageBaseType(new MatrixType(rng_spc.mapper().size(), src_spc.mapper().size(), pattern(rng_spc, src_spc, grd_vw)))
    , AssemblerBaseType(rng_spc, src_spc, grd_vw)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , coupling_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
    , boundary_operator_(over_integrate, diffusion_factor_, diffusion_tensor_)
    , coupling_assembler_(coupling_operator_)
    , boundary_assembler_(boundary_operator_)
    , assembled_(false)
  {
    setup();
  }

  const GridViewType& grid_view() const
  {
    return AssemblerBaseType::grid_view();
  }

  const RangeSpaceType& range_space() const
  {
    return AssemblerBaseType::test_space();
  }

  const SourceSpaceType& source_space() const
  {
    return AssemblerBaseType::ansatz_space();
  }

  MatrixType& matrix()
  {
    return StorageBaseType::storage_access();
  }

  const MatrixType& matrix() const
  {
    return StorageBaseType::storage_access();
  }

  void assemble()
  {
    if (!assembled_) {
      AssemblerBaseType::assemble();
      assembled_ = true;
    }
  } // ... assemble()

private:
  void setup()
  {
    this->add(coupling_assembler_, matrix(), new DSG::ApplyOn::InnerIntersectionsPrimally< GridViewType >());
    this->add(boundary_assembler_, matrix(), new DSG::ApplyOn::BoundaryIntersections< GridViewType >());
  }

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const CouplingOperatorType coupling_operator_;
  const BoundaryOperatorType boundary_operator_;
  const CouplingAssemblerType coupling_assembler_;
  const BoundaryAssemblerType boundary_assembler_;
  bool assembled_;
}; // class EllipticSWIPDGPenaltyAssemblable


} // namespace Products
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_PRODUCTS_BASE_HH
