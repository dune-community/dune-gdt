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
#include <dune/gdt/localevaluation/elliptic-ipdg.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/codim1.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/codim1.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forward
template <class DiffusionFactorType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp = typename SourceSpaceImp::GridViewType, class DiffusionTensorType = void>
class EllipticSWIPDG;


namespace internal {


template <class DiffusionFactorType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp,
          class DiffusionTensorType = void>
class EllipticSWIPDGTraits
{
  static_assert(Stuff::is_localizable_function<DiffusionFactorType>::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function<DiffusionTensorType>::value,
                "DiffusionTensorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::LA::is_matrix<MatrixImp>::value, "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(is_space<SourceSpaceImp>::value, "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(is_space<RangeSpaceImp>::value, "RangeSpaceImp has to be derived from SpaceInterface!");

public:
  typedef EllipticSWIPDG<DiffusionFactorType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp,
                         DiffusionTensorType> derived_type;
  typedef MatrixImp MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef GridViewImp GridViewType;
}; // class EllipticSWIPDGTraits


template <class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void>
{
  static_assert(Stuff::is_localizable_function<DiffusionType>::value,
                "DiffusionType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::LA::is_matrix<MatrixImp>::value, "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(is_space<SourceSpaceImp>::value, "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(is_space<RangeSpaceImp>::value, "RangeSpaceImp has to be derived from SpaceInterface!");

public:
  typedef EllipticSWIPDG<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void> derived_type;
  typedef MatrixImp MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef GridViewImp GridViewType;
}; // class EllipticSWIPDGTraits


} // namespace internal


template <class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class EllipticSWIPDG<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void>
    : Stuff::Common::StorageProvider<MatrixImp>,
      public Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp,
                                                                   RangeSpaceImp, GridViewImp, void>>,
      public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef Stuff::Common::StorageProvider<MatrixImp> StorageProvider;
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> AssemblerBaseType;
  typedef Operators::MatrixBased<internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp,
                                                                GridViewImp, void>> OperatorBaseType;

  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> VolumeOperatorType;
  typedef LocalAssembler::Codim0Matrix<VolumeOperatorType> VolumeAssemblerType;
  typedef LocalOperator::Codim1CouplingIntegral<LocalEvaluation::SWIPDG::Inner<DiffusionType>> CouplingOperatorType;
  typedef LocalAssembler::Codim1CouplingMatrix<CouplingOperatorType> CouplingAssemblerType;
  typedef LocalOperator::Codim1BoundaryIntegral<LocalEvaluation::SWIPDG::BoundaryLHS<DiffusionType>>
      DirichletBoundaryOperatorType;
  typedef LocalAssembler::Codim1BoundaryMatrix<DirichletBoundaryOperatorType> DirichletBoundaryAssemblerType;
  typedef typename MatrixImp::ScalarType ScalarType;

public:
  typedef internal::EllipticSWIPDGTraits<DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void>
      Traits;
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
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(matrix)
    , OperatorBaseType(this->storage_access(), source_space, range_space, grid_view)
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
    setup();
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info,
                 const SourceSpaceType& source_space, const RangeSpaceType& range_space, const GridViewType& grid_view,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(new MatrixType(range_space.mapper().size(), source_space.mapper().size(),
                                     pattern(range_space, source_space, grid_view)))
    , OperatorBaseType(this->storage_access(), source_space, range_space, grid_view)
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
    setup();
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info, MatrixType& matrix,
                 const SourceSpaceType& source_space, const RangeSpaceType& range_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(matrix)
    , OperatorBaseType(this->storage_access(), source_space, range_space)
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
    setup();
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info,
                 const SourceSpaceType& source_space, const RangeSpaceType& range_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(
          new MatrixType(range_space.mapper().size(), source_space.mapper().size(), pattern(range_space, source_space)))
    , OperatorBaseType(this->storage_access(), source_space, range_space)
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
    setup();
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info, MatrixType& matrix,
                 const SourceSpaceType& source_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(matrix)
    , OperatorBaseType(this->storage_access(), source_space)
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
    setup();
  }

  EllipticSWIPDG(const DiffusionType& diffusion, const BoundaryInfoType& boundary_info,
                 const SourceSpaceType& source_space,
                 const ScalarType beta = LocalEvaluation::SIPDG::internal::default_beta(GridViewType::dimension))
    : StorageProvider(new MatrixType(source_space.mapper().size(), source_space.mapper().size(), pattern(source_space)))
    , OperatorBaseType(this->storage_access(), source_space)
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

  const DiffusionType& diffusion_;
  const BoundaryInfoType& boundary_info_;
  const VolumeOperatorType volume_operator_;
  const VolumeAssemblerType volume_assembler_;
  const CouplingOperatorType coupling_operator_;
  const CouplingAssemblerType coupling_assembler_;
  const DirichletBoundaryOperatorType dirichlet_boundary_operator_;
  const DirichletBoundaryAssemblerType dirichlet_boundary_assembler_;
}; // class EllipticSWIPDG


/// \todo use matrix as first template parameter, dro /*matrix*/
/// \todo return by value, implement move ctor
template <class DF, class M, class S>
std::unique_ptr<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType>>
make_elliptic_swipdg(const DF& diffusion_factor,
                     const Stuff::Grid::BoundaryInfoInterface<typename S::GridViewType::Intersection>& boundary_info,
                     const M& /*matrix*/, const S& space)
{
  return Stuff::Common::make_unique<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType>>(
      diffusion_factor, boundary_info, space);
} // ... make_elliptic_swipdg(...)

template <class DF, class M, class S>
std::unique_ptr<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType>>
make_elliptic_swipdg(const DF& diffusion_factor,
                     const Stuff::Grid::BoundaryInfoInterface<typename S::GridViewType::Intersection>& boundary_info,
                     M& matrix, const S& space)
{
  return Stuff::Common::make_unique<EllipticSWIPDG<DF, M, S, S, typename S::GridViewType>>(
      diffusion_factor, boundary_info, matrix, space);
} // ... make_elliptic_swipdg(...)


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_SWIPDG_HH
