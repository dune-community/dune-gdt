// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_CG_HH
#define DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_CG_HH

#include <dune/gdt/playground/localevaluation/elliptic.hh>
#include <dune/gdt/operators/elliptic-cg.hh>

namespace Dune {
namespace GDT {
namespace Operators {


template <class DiffusionFactorType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp,
          class DiffusionTensorType>
class EllipticCG
    : public Operators::MatrixBased<internal::EllipticCGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp,
                                                               RangeSpaceImp, GridViewImp, DiffusionTensorType>>,
      public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> AssemblerBaseType;
  typedef Operators::MatrixBased<internal::EllipticCGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp,
                                                            RangeSpaceImp, GridViewImp, DiffusionTensorType>>
      OperatorBaseType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionFactorType, DiffusionTensorType>>
      LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix<LocalOperatorType> LocalAssemblerType;

public:
  typedef internal::EllipticCGTraits<DiffusionFactorType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp,
                                     DiffusionTensorType> Traits;

  typedef typename Traits::MatrixType MatrixType;
  typedef typename Traits::SourceSpaceType SourceSpaceType;
  typedef typename Traits::RangeSpaceType RangeSpaceType;
  typedef typename Traits::GridViewType GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space, const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  EllipticCG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
             MatrixType& matrix, const SourceSpaceType& source_space, const RangeSpaceType& range_space,
             const GridViewType& grid_view)
    : OperatorBaseType(matrix, source_space, range_space, grid_view)
    , AssemblerBaseType(range_space, grid_view, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_operator_(diffusion_factor_, diffusion_tensor_)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
             MatrixType& matrix, const SourceSpaceType& source_space, const RangeSpaceType& range_space)
    : OperatorBaseType(matrix, source_space, range_space)
    , AssemblerBaseType(range_space, source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_operator_(diffusion_factor_, diffusion_tensor_)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionFactorType& diffusion_factor, const DiffusionTensorType& diffusion_tensor,
             MatrixType& matrix, const SourceSpaceType& source_space)
    : OperatorBaseType(matrix, source_space)
    , AssemblerBaseType(source_space)
    , diffusion_factor_(diffusion_factor)
    , diffusion_tensor_(diffusion_tensor)
    , local_operator_(diffusion_factor_, diffusion_tensor_)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  virtual ~EllipticCG()
  {
  }

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  const DiffusionFactorType& diffusion_factor_;
  const DiffusionTensorType& diffusion_tensor_;
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
}; // class EllipticCG


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_PLAYGROUND_OPERATORS_ELLIPTIC_CG_HH
