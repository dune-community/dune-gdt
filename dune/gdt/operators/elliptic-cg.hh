// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_ELLIPTIC_CG_HH
#define DUNE_GDT_OPERATORS_ELLIPTIC_CG_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune {
namespace GDT {
namespace Operators {


// forwards
template< class DiffusionFactorType
        , class MatrixImp
        , class SourceSpaceImp
        , class RangeSpaceImp = SourceSpaceImp
        , class GridViewImp = typename SourceSpaceImp::GridViewType
        , class DiffusionTensorType = void >
class EllipticCG;


namespace internal {


template< class DiffusionFactorType
        , class MatrixImp
        , class SourceSpaceImp
        , class RangeSpaceImp
        , class GridViewImp
        , class DiffusionTensorType = void >
class EllipticCGTraits
{
  static_assert(Stuff::is_localizable_function< DiffusionFactorType >::value,
                "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_localizable_function< DiffusionTensorType >::value
                || std::is_same< void, DiffusionTensorType>::value,
                "DiffusionTensorType has to be void or derived from Stuff::LocalizableFunctionInterface!");
  static_assert(Stuff::is_matrix< MatrixImp >::value,
                "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
  static_assert(is_space< SourceSpaceImp >::value, "SourceSpaceImp has to be derived from SpaceInterface!");
  static_assert(is_space< RangeSpaceImp >::value,  "RangeSpaceImp has to be derived from SpaceInterface!");
public:
  typedef EllipticCG< DiffusionFactorType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, DiffusionTensorType >
      derived_type;
  typedef MatrixImp       MatrixType;
  typedef SourceSpaceImp  SourceSpaceType;
  typedef RangeSpaceImp   RangeSpaceType;
  typedef GridViewImp     GridViewType;
}; // class EllipticCGTraits


} // namespace internal


template< class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp >
class EllipticCG< DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void >
  : Stuff::Common::StorageProvider< MatrixImp >
  , public Operators::MatrixBased< internal::EllipticCGTraits< DiffusionType, MatrixImp
                                                             , SourceSpaceImp, RangeSpaceImp, GridViewImp, void > >
  , public SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >
{
  typedef Stuff::Common::StorageProvider< MatrixImp > StorageProvider;
  typedef SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp > AssemblerBaseType;
  typedef Operators::MatrixBased< internal::EllipticCGTraits< DiffusionType, MatrixImp, SourceSpaceImp
                                                            , RangeSpaceImp, GridViewImp > > OperatorBaseType;
  typedef LocalOperator::Codim0Integral< LocalEvaluation::Elliptic< DiffusionType > > LocalOperatorType;
  typedef LocalAssembler::Codim0Matrix< LocalOperatorType >                           LocalAssemblerType;
public:
  typedef internal::EllipticCGTraits< DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void >
      Traits;

  typedef typename Traits::MatrixType       MatrixType;
  typedef typename Traits::SourceSpaceType  SourceSpaceType;
  typedef typename Traits::RangeSpaceType   RangeSpaceType;
  typedef typename Traits::GridViewType     GridViewType;

  using OperatorBaseType::pattern;

  static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                   const SourceSpaceType& source_space,
                                                   const GridViewType& grid_view)
  {
    return range_space.compute_volume_pattern(grid_view, source_space);
  }

  EllipticCG(const DiffusionType& diffusion,
             MatrixType& mtrx,
             const SourceSpaceType& src_spc,
             const RangeSpaceType& rng_spc,
             const GridViewType& grid_view)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_spc, rng_spc, grid_view)
    , AssemblerBaseType(rng_spc, grid_view, src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion,
             const SourceSpaceType& src_spc,
             const RangeSpaceType& rng_spc,
             const GridViewType& grid_view)
    : StorageProvider(new MatrixType(rng_spc.mapper().size(),
                                     src_spc.mapper().size(),
                                     pattern(rng_spc, src_spc, grid_view)))
    , OperatorBaseType(this->storage_access(), src_spc, rng_spc, grid_view)
    , AssemblerBaseType(rng_spc, grid_view, src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion,
             MatrixType& mtrx,
             const SourceSpaceType& src_spc,
             const RangeSpaceType& rng_spc)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_spc, rng_spc)
    , AssemblerBaseType(rng_spc, src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion,
             const SourceSpaceType& src_spc,
             const RangeSpaceType& rng_spc)
    : StorageProvider(new MatrixType(rng_spc.mapper().size(),
                                     src_spc.mapper().size(),
                                     pattern(rng_spc, src_spc)))
    , OperatorBaseType(this->storage_access(), src_spc, rng_spc)
    , AssemblerBaseType(rng_spc, src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion,
             MatrixType& mtrx,
             const SourceSpaceType& src_spc)
    : StorageProvider(mtrx)
    , OperatorBaseType(this->storage_access(), src_spc)
    , AssemblerBaseType(src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  EllipticCG(const DiffusionType& diffusion,
             const SourceSpaceType& src_spc)
    : StorageProvider(new MatrixType(src_spc.mapper().size(), src_spc.mapper().size(), pattern(src_spc)))
    , OperatorBaseType(this->storage_access(), src_spc)
    , AssemblerBaseType(src_spc)
    , diffusion_(diffusion)
    , local_operator_(diffusion_)
    , local_assembler_(local_operator_)
    , assembled_(false)
  {
    this->add(local_assembler_, this->matrix());
  }

  virtual ~EllipticCG() {}

  virtual void assemble() override final
  {
    if (!assembled_) {
      AssemblerBaseType::assemble(true);
      assembled_ = true;
    }
  } // ... assemble(...)

private:
  const DiffusionType& diffusion_;
  const LocalOperatorType local_operator_;
  const LocalAssemblerType local_assembler_;
  bool assembled_;
}; // class EllipticCG


} // namespace Operators
} // namespace GDT
} // namespace Dune

#endif // DUNE_GDT_OPERATORS_ELLIPTIC_CG_HH
