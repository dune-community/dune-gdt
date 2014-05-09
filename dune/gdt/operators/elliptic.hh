// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_HH

#include <type_traits>

#include <dune/stuff/la/container/interfaces.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/common/memory.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/assembler/local/codim0.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forward, to be used in the traits
template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp                                                              = typename SourceSpaceImp::GridViewType>
class EllipticCG;


template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp = SourceSpaceImp,
          class GridViewImp                                                              = typename SourceSpaceImp::GridViewType>
class EllipticCGTraits
{
  static_assert(std::is_base_of<Stuff::LocalizableFunctionInterface<
                                    typename DiffusionImp::EntityType, typename DiffusionImp::DomainFieldType,
                                    DiffusionImp::dimDomain, typename DiffusionImp::RangeFieldType,
                                    DiffusionImp::dimRange, DiffusionImp::dimRangeCols>,
                                DiffusionImp>::value,
                "DiffusionImp has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(std::is_base_of<Stuff::LA::MatrixInterface<typename MatrixImp::Traits>, MatrixImp>::value,
                "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename SourceSpaceImp::Traits>, SourceSpaceImp>::value,
                "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(std::is_base_of<SpaceInterface<typename RangeSpaceImp::Traits>, RangeSpaceImp>::value,
                "RangeSpaceImp has to be derived from SpaceInterface!");

public:
  typedef EllipticCG<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> derived_type;
  typedef MatrixImp MatrixType;
  typedef SourceSpaceImp SourceSpaceType;
  typedef RangeSpaceImp RangeSpaceType;
  typedef GridViewImp GridViewType;

private:
  typedef DiffusionImp DiffusionType;

public:
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> LocalOperatorType;

private:
  friend class EllipticCG<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp>;
}; // class EllipticTraits


template <class DiffusionImp, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp>
class EllipticCG : public Operators::MatrixBasedBase<internal::EllipticCGTraits<DiffusionImp, MatrixImp, SourceSpaceImp,
                                                                                RangeSpaceImp, GridViewImp>>,
                   public SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp>
{
  typedef SystemAssembler<RangeSpaceImp, GridViewImp, SourceSpaceImp> AssemblerBaseType;
  typedef Operators::MatrixBasedBase<internal::EllipticCGTraits<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp,
                                                                GridViewImp>> OperatorBaseType;

  typedef DiffusionImp DiffusionType;
  typedef LocalOperator::Codim0Integral<LocalEvaluation::Elliptic<DiffusionType>> LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix<LocalOperatorType> LocalAssemblerType;

public:
  typedef internal::EllipticCGTraits<DiffusionImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp> Traits;

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

  EllipticCG(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space,
             const RangeSpaceType& range_space, const GridViewType& grid_view)
    : OperatorBaseType(matrix, source_space, range_space, grid_view)
    , AssemblerBaseType(range_space, grid_view, source_space)
    , local_operator_(diffusion)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space,
             const RangeSpaceType& range_space)
    : OperatorBaseType(matrix, source_space, range_space)
    , AssemblerBaseType(range_space, source_space)
    , local_operator_(diffusion)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion, MatrixType& matrix, const SourceSpaceType& source_space)
    : OperatorBaseType(matrix, source_space)
    , AssemblerBaseType(source_space)
    , local_operator_(diffusion)
    , local_assembler_(local_operator_)
  {
    this->add(local_assembler_, this->matrix());
  }

  virtual void assemble() DS_OVERRIDE DS_FINAL
  {
    AssemblerBaseType::assemble();
  }

private:
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
}; // class EllipticCG


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_HH
