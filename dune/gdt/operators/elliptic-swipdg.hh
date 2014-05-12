// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_SWIPDG_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_SWIPDG_HH

#include <type_traits>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>
#include <dune/stuff/grid/boundaryinfo.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/swipdg.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forwards
template <class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp                                                               = typename SourceSpaceImp::GridViewType>
class EllipticSWIPDG;


namespace internal {


template <class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class EllipticSWIPDGTraits
{
  static_assert(std::is_base_of<Stuff::LocalizableFunctionInterface<
                                    typename DiffusionType::EntityType, typename DiffusionType::DomainFieldType,
                                    DiffusionType::dimDomain, typename DiffusionType::RangeFieldType,
                                    DiffusionType::dimRange, DiffusionType::dimRangeCols>,
                                DiffusionType>::value,
                "DiffusionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixImp::Traits>, MatrixImp>::value,
                "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceImp::Traits>, SourceSpaceImp>::value,
                "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceImp::Traits>, RangeSpaceImp>::value,
                "RangeSpaceImp has to be derived from SpaceInterface!");

public:
  typedef EllipticSWIPDG<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> derived_type;
  typedef MatrixImp MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef GridViewImp GridViewType;
}; // class EllipticSWIPDGTraits


} // namespace internal


template <class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class EllipticSWIPDG
    : public Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp,
                                                                   RangeSpaceImp, GridViewImp>>,
      public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> AssemblerBaseType;
  typedef Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp,
                                                                GridViewImp>> OperatorBaseType;

  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> VolumeOperatorType;
  typedef LocalAssembler::Codim0Matrix<VolumeOperatorType> VolumeAssemblerType;
  typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDG::Inner<DiffusionType>> CouplingOperatorType;
  typedef LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> CouplingAssemblerType;

  typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SWIPDG::BoundaryLHS<DiffusionType>>
      DirichletBoundaryOperatorType;
  typedef LocalAssembler::Codim1BoundaryMatrix<DirichletBoundaryOperatorType> DirichletBoundaryAssemblerType;

  typedef typename MatrixImp::ScalarType ScalarType;

public:
  typedef internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> Traits;

  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;

  typedef Stuff::Grid::BoundaryInfoInterface<typename GridViewType::Intersection> BoundaryInfoType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_face_and_volume_pattern(grid_view, source_space);
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info, MatrixType& matrix,
                 const SourceSpaceType& source_space, const RangeSpaceType& range_space, const GridViewType& grid_view,
                 const ScalarType beta = 1.0)
    : OperatorBaseType(matrix, source_space, range_space, grid_view)
    , AssemblerBaseType(range_space, grid_view, source_space)
    , diffusion_(diffusion)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    this->add(volume_assembler_, this->matrix());
    this->add(coupling_assembler_, this->matrix(), new ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(dirichlet_boundary_assembler_,
              this->matrix(),
              new ApplyOn::DirichletIntersections<GridViewType>(boundary_info_));
  } // EllipticSWIPDG(...)

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info, MatrixType& matrix,
                 const SourceSpaceType& source_space, const RangeSpaceType& range_space, const ScalarType beta = 1.0)
    : OperatorBaseType(matrix, source_space, range_space)
    , AssemblerBaseType(range_space, source_space)
    , diffusion_(diffusion)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    this->add(volume_assembler_, this->matrix());
    this->add(coupling_assembler_, this->matrix(), new ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(dirichlet_boundary_assembler_,
              this->matrix(),
              new ApplyOn::DirichletIntersections<GridViewType>(boundary_info_));
  } // EllipticSWIPDG(...)

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info, MatrixType& matrix,
                 const SourceSpaceType& source_space, const ScalarType beta = 1.0)
    : OperatorBaseType(matrix, source_space)
    , AssemblerBaseType(source_space)
    , diffusion_(diffusion)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    this->add(volume_assembler_, this->matrix());
    this->add(coupling_assembler_, this->matrix(), new ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(dirichlet_boundary_assembler_,
              this->matrix(),
              new ApplyOn::DirichletIntersections<GridViewType>(boundary_info_));
  } // EllipticSWIPDG(...)

  virtual ~EllipticSWIPDG()
  {
  }

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  const DiffusionType& diffusion_;
  const BoundaryInfoType& boundary_info_;
  const VolumeOperatorType volume_operator_;
  const VolumeAssemblerType volume_assembler_;
  const CouplingOperatorType coupling_operator_;
  const CouplingAssemblerType coupling_assembler_;
  const DirichletBoundaryOperatorType dirichlet_boundary_operator_;
  const DirichletBoundaryAssemblerType dirichlet_boundary_assembler_;
}; // class EllipticSWIPDG


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_SWIPDG_HH
