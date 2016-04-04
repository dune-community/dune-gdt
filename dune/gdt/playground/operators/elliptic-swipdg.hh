// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_SWIPDG_HH
#define DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_SWIPDG_HH

#include <dune/stuff/common/memory.hh>

#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/operators/elliptic-swipdg.hh>

namespace Dune {
namespace GDT {
namespace Operators {


template <class DiffusionFactorType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp,
          class DiffusionTensorType>
class EllipticSWIPDG
    : Stuff::Common::StorageProvider<MatrixImp>,
      public Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp,
                                                                   RangeSpaceImp, GridViewImp, DiffusionTensorType>>,
      public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef Stuff::Common::StorageProvider<MatrixImp> StorageBaseType;
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> AssemblerBaseType;
  typedef Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp,
                                                                RangeSpaceImp, GridViewImp, DiffusionTensorType>>
      OperatorBaseType;

  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionFactorType, DiffusionTensorType>>
      VolumeOperatorType;
  typedef LocalAssembler::Codim0Matrix<VolumeOperatorType> VolumeAssemblerType;
  typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDG::Inner<DiffusionFactorType,
                                                                               DiffusionTensorType>>
      CouplingOperatorType;
  typedef LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> CouplingAssemblerType;
  typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SWIPDG::BoundaryLHS<DiffusionFactorType,
                                                                                     DiffusionTensorType>>
      DirichletBoundaryOperatorType;
  typedef LocalAssembler::Codim1BoundaryMatrix<DirichletBoundaryOperatorType> DirichletBoundaryAssemblerType;

  typedef typename MatrixImp::ScalarType ScalarType;

public:
  typedef internal::EllipticSWIPDGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp,
                                         DiffusionTensorType> Traits;

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

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, MatrixType& matrix, const SourceSpaceType& source_space,
                 const RangeSpaceType& range_space, const GridViewType& grid_view,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(matrix)
    , OperatorBaseType(this->storage_access(), source_space, range_space, grid_view)
    , AssemblerBaseType(range_space, grid_view, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, const SourceSpaceType& source_space,
                 const RangeSpaceType& range_space, const GridViewType& grid_view,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(new MatrixType(range_space.mapper().size(), source_space.mapper().size(),
                                     pattern(range_space, source_space, grid_view)))
    , OperatorBaseType(this->storage_access(), source_space, range_space, grid_view)
    , AssemblerBaseType(range_space, grid_view, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, MatrixType& matrix, const SourceSpaceType& source_space,
                 const RangeSpaceType& range_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(matrix)
    , OperatorBaseType(this->storage_access(), source_space, range_space)
    , AssemblerBaseType(range_space, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, const SourceSpaceType& source_space,
                 const RangeSpaceType& range_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(
          new MatrixType(range_space.mapper().size(), source_space.mapper().size(), pattern(range_space, source_space)))
    , OperatorBaseType(this->storage_access(), source_space, range_space)
    , AssemblerBaseType(range_space, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, MatrixType& matrix, const SourceSpaceType& source_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(matrix)
    , OperatorBaseType(this->storage_access(), source_space)
    , AssemblerBaseType(source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  EllipticSWIPDG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
                 const BoundaryInfoType& boundary_info, const SourceSpaceType& source_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageBaseType(new MatrixType(source_space.mapper().size(), source_space.mapper().size(), pattern(source_space)))
    , OperatorBaseType(this->storage_access(), source_space)
    , AssemblerBaseType(source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , boundary_info_(boundary_info)
    , volume_operator_(diffusion_factor_, diffusion_tensor_)
    , volume_assembler_(volume_operator_)
    , coupling_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , coupling_assembler_(coupling_operator_)
    , dirichlet_boundary_operator_(diffusion_factor_, diffusion_tensor_, beta)
    , dirichlet_boundary_assembler_(dirichlet_boundary_operator_)
  {
    setup();
  }

  virtual ~EllipticSWIPDG()
  {
  }

  virtual void assemble() override final
  {
    AssemblerBaseType::assemble();
  }

private:
  void setup()
  {
    this->add(volume_assembler_, this->matrix());
    this->add(
        coupling_assembler_, this->matrix(), new Stuff::Grid::ApplyOn::InnerIntersectionsPrimally<GridViewType>());
    this->add(dirichlet_boundary_assembler_,
              this->matrix(),
              new Stuff::Grid::ApplyOn::DirichletIntersections<GridViewType>(boundary_info_));
  } // ... setup(...)

  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const BoundaryInfoType& boundary_info_;
  const VolumeOperatorType volume_operator_;
  const VolumeAssemblerType volume_assembler_;
  const CouplingOperatorType coupling_operator_;
  const CouplingAssemblerType coupling_assembler_;
  const DirichletBoundaryOperatorType dirichlet_boundary_operator_;
  const DirichletBoundaryAssemblerType dirichlet_boundary_assembler_;
}; // class EllipticSWIPDG


template <class DF, class DT, class M, class S>
std::unique_ptr<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType, DT>>
make_elliptic_swipdg(const DF& diffusion_factor, const DT& diffusion_tensor,
                     const Stuff::Grid::BoundaryInfoInterface<typename S::GridViewType::Intersection>& boundary_info,
                     const M& /*matrix*/, const S& space)
{
  return Stuff::Common::make_unique<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType, DT>>(
      diffusion_factor, diffusion_tensor, boundary_info, space);
} // ... make_elliptic_swipdg(...)

template <class DF, class DT, class M, class S>
std::unique_ptr<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType, DT>>
make_elliptic_swipdg(const DF& diffusion_factor, const DT& diffusion_tensor,
                     const Stuff::Grid::BoundaryInfoInterface<typename S::GridViewType::Intersection>& boundary_info,
                     M& matrix, const S& space)
{
  return Stuff::Common::make_unique<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType, DT>>(
      diffusion_factor, diffusion_tensor, boundary_info, matrix, space);
} // ... make_elliptic_swipdg(...)


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_SWIPDG_HH
