// This file is part of the dune-gdt project:
//   http://users.dune-project.org/projects/dune-gdt
// Copyright holders: Felix Schindler
// License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)

#ifndef DUNE_GDT_OPERATORS_VIRTUAL_REFINED_HH
#define DUNE_GDT_OPERATORS_VIRTUAL_REFINED_HH

#include <type_traits>

#include <dune/stuff/common/memory.hh>
#include <dune/stuff/functions/interfaces.hh>
#include <dune/stuff/la/container/interfaces.hh>

#include <dune/gdt/spaces/interface.hh>
#include <dune/gdt/localevaluation/elliptic.hh>
#include <dune/gdt/localoperator/codim0.hh>
#include <dune/gdt/localoperator/virtual_codim0.hh>
#include <dune/gdt/assembler/local/codim0.hh>
#include <dune/gdt/assembler/local/virtual_codim0.hh>
#include <dune/gdt/assembler/system.hh>

#include "base.hh"

namespace Dune
{
namespace GDT
{
namespace Operators
{


// forward
template< class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp,
          class GridViewImp, class DiffusionTensorType = void>
class VirtualRefinedEllipticCG;


namespace internal
{


  template< class DiffusionFactorImp
          , class MatrixImp
          , class SourceSpaceImp
          , class RangeSpaceImp
          , class GridViewImp
          , class DiffusionTensorType = void >
  class VirtualRefinedEllipticCGTraits
  {
    static_assert(Stuff::is_localizable_function< DiffusionFactorImp >::value,
                  "DiffusionFactorType has to be derived from Stuff::LocalizableFunctionInterface!");
    static_assert(Stuff::is_localizable_function< DiffusionTensorType >::value
                  || std::is_same< void, DiffusionTensorType>::value,
                  "DiffusionTensorType has to be void or derived from Stuff::LocalizableFunctionInterface!");
    static_assert(Stuff::LA::is_matrix< MatrixImp >::value,
                  "MatrixImp has to be derived from Stuff::LA::MatrixInterface!");
    static_assert(is_space< SourceSpaceImp >::value, "SourceSpaceImp has to be derived from SpaceInterface!");
    static_assert(is_space< RangeSpaceImp >::value,  "RangeSpaceImp has to be derived from SpaceInterface!");
  public:
    typedef VirtualRefinedEllipticCG< DiffusionFactorImp, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, DiffusionTensorType >
        derived_type;
    typedef MatrixImp       MatrixType;
    typedef SourceSpaceImp  SourceSpaceType;
    typedef RangeSpaceImp   RangeSpaceType;
    typedef GridViewImp     GridViewType;
    typedef DiffusionFactorImp DiffusionFactorType;
  }; // class EllipticCGTraits


  } // namespace internal


  template< class DiffusionType, class MatrixImp, class SourceSpaceImp, class RangeSpaceImp, class GridViewImp >
  class VirtualRefinedEllipticCG< DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void >
    : Stuff::Common::StorageProvider< MatrixImp >
    , public Operators::MatrixBased< internal::VirtualRefinedEllipticCGTraits< DiffusionType, MatrixImp
                                                               , SourceSpaceImp, RangeSpaceImp, GridViewImp, void > >
    , public SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >
  {
    typedef Stuff::Common::StorageProvider< MatrixImp >                                        StorageProvider;
    typedef SystemAssembler< RangeSpaceImp, GridViewImp, SourceSpaceImp >                      AssemblerBaseType;
    typedef Operators::MatrixBased< internal::VirtualRefinedEllipticCGTraits< DiffusionType, MatrixImp, SourceSpaceImp
                                                              , RangeSpaceImp, GridViewImp > > OperatorBaseType;
    typedef LocalOperator::VirtualRefinedCodim0Integral< LocalEvaluation::Elliptic< DiffusionType > >        VirtualRefinedLocalOperatorType;
    typedef LocalAssembler::Codim0Matrix< VirtualRefinedLocalOperatorType >                                  LocalAssemblerType;
  public:
    typedef internal::VirtualRefinedEllipticCGTraits< DiffusionType, MatrixImp, SourceSpaceImp, RangeSpaceImp, GridViewImp, void >
        Traits;

    typedef typename Traits::MatrixType      MatrixType;
    typedef typename Traits::SourceSpaceType SourceSpaceType;
    typedef typename Traits::RangeSpaceType  RangeSpaceType;
    typedef typename Traits::GridViewType    GridViewType;

    using OperatorBaseType::pattern;

    static Stuff::LA::SparsityPatternDefault pattern(const RangeSpaceType& range_space,
                                                     const SourceSpaceType& source_space,
                                                     const GridViewType& grid_view)
    {
      return range_space.compute_volume_pattern(grid_view, source_space);
    }

  VirtualRefinedEllipticCG(const DiffusionType& diffusion, MatrixType& mtrx, const SourceSpaceType& src_spc)
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
  virtual void assemble() override final
  {
    if (!assembled_) {
      AssemblerBaseType::assemble(true);
      assembled_ = true;
    }
  } // ... assemble(...)
private:
  const DiffusionType& diffusion_;
  const VirtualRefinedLocalOperatorType local_operator_;
  LocalAssemblerType local_assembler_;
  bool assembled_;
}; // class VirtualRefinedEllipticCG

} // namespace Operators
} // namespace GDT
} // namespace Dune


#endif // VIRTUAL_REFINED_HH
